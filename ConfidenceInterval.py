#!/usr/bin/env python

from __future__ import division

import inspect
import logging
logging.basicConfig()

import random
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import poisson
from scipy.stats import norm

from scipy import integrate
from scipy import optimize

from math import log

class variable(object):
    """ A class to store a name, value, and range

    """

    def __init__(self, name, var_min=-np.inf, var_max=np.inf, num_points=1000):
        self.name = name
        self.min = var_min
        self.max = var_max
        self.num_points = num_points
        #self.val = 0.0

    def linspace(self):
        return np.linspace(self.min, self.max, self.num_points)

    pass


class likelihood(object):
    """ A function which evaluates a pdf on data

    A likelihood class takes a pdf: pdf(data, *params)
    One can specify the data argument and the various
    parameters by name, or by supplying 'variables'
    (which then give those arguments min, max and num_points)
    
    """

    # Create an internal logger
    logging = logging.getLogger("likelihood")

    def __init__(self, pdf, data=None, params=None):

        # Set the pdf
        self.pdf = pdf
        self.args = {}

        # Use inspection to find all arguments of the pdf
        func_spec = inspect.getargspec(self.pdf)
        (all_arguments, all_defaults) = (func_spec.args, func_spec.defaults)

        # Get the defaults for any arguments that have them
        default_dict = {}
        for name, val in zip( reversed(all_arguments), reversed(all_defaults) ):
            default_dict[name] = val

        # Determine the data var and save an internal 'variable'
        # If none is supplied, assume it's the 0th argument
        # If a string is supplied, create a new variable
        if data == None:
            data = all_arguments[0]
        if isinstance(data, variable):
            self.data = data
        else:
            self.data = variable(data)
        if self.data.name not in all_arguments:
            print "Error: Supplied data var is not an argument of the supplied pdf",
            print all_arguments
            raise Exception("InvalidDataVar")

        # And create an attribute for easy access
        if data.name in default_dict:
            setattr(self, data.name, default_dict[data.name])
        else:
            setattr(self, data.name, 0.0)
        #self.param_dict[data.name] = data
            
        # Get the list of parameters,
        # create variables based on them,
        # and store that list of variables
        self.param_dict = {}
        param_list = [arg for arg in func_spec.args if arg != self.data.name]
        for param in param_list:

            # Check if the parameter matches one in the
            # supplied list of parameter 'variables'
            matching_var_list = [var for var in params if var.name==param]
            
            param_var = None

            if len(matching_var_list)==0:
                self.logging.debug("Creating new param with name: " + param)
                param_var = variable(param)
                # If the parameter has a default in the function definition,
                # Set that default here
                # MUST SET: param_var.val = param.defaults
            elif len(matching_var_list)==1:
                param_var = matching_var_list[0]
            else:
                print "Error: More than one parameter variable supplied ",
                print "with the name: ", param
                raise Exception("ParamVariable")
            self.param_dict[param_var.name] = param_var

            # And create an attribute for easy access
            if param_var.name in default_dict:
                setattr(self, param_var.name, default_dict[param_var.name])
            else:
                setattr(self, param_var.name, 0.0)
            setattr(self, param_var.name + "_var", param_var)

        self.args.update(self.param_dict)
        self.args[self.data.name] = self.data

        self.norm = 1.0
        self.normalization_cache = {}
        self.minimization_cache = {}
        self.nll_cache = {}


    def get_data(self):
        """ Get the current value of the data
        """
        data_name = self.data.name
        return getattr(self, data_name)


    def set_data(self, val):
        """ Get the current value of the data
        """
        data_name = self.data.name
        return setattr(self, data_name, val)


    def param_state(self):
        """ Return a dict with the current state only the parameters

        """
        current_state = {}
        #current_state[self.data.name] = self.data.val
        for name, param in self.param_dict.iteritems():
            current_state[param.name] = getattr(self, param.name)
        return current_state


    def total_state(self):
        """ Return a dict with the current state including data and all parameters

        """
        current_state = {}
        #current_state[self.data.name] = self.data.val
        for name, param in self.param_dict.iteritems():
            current_state[param.name] = getattr(self, param.name)
        current_state[self.data.name] = self.get_data()
        return current_state


    def set_state(self, **kwargs):
        """ Set the state based on the values of the arguments

        Return any args that aren't parameters of the likelihood
        """

        for (arg, val) in kwargs.iteritems():
            if hasattr(self, arg):
                setattr(self, arg, val)
            else:
                print "Error: Cannot set state argument: ", arg,
                print " in likelihood, it does not exist"
                raise Exception("SetState")
            pass
        return


    def _eval_raw(self, **kwargs):
        """ Get the current state of the likelihood
        without any normalization

        Any values can be set via kwargs, and the evaluation
        will take place after those kwargs are set
        """

        self.set_state(**kwargs)
        current_state = self.total_state()
        #return self.pdf(self.get_data(), **current_state)
        pdf_val = self.pdf(**current_state)
        if pdf_val < 0:
            print "Error: Pdf is 0 at state: ", current_state
            raise Exception("PdfVal")
        return pdf_val
        

    def eval(self, **kwargs):
        """ Val of pdf based on the current state,
        Evaluated on the given data point
        This includes normalization, which is cached

        The current state can be set via kwargs, and
        the normalization will take place after the
        current state is set.

        """
        self.set_state(**kwargs)

        # We have to normalize, but we should be sure
        # to restore the state after, so we aren't effected
        # by the random data points used to evaluate the integral
        #current_state = self.total_state()
        #current_data = self.get_data()
        self.normalize()
        #self.set_data(current_data)
        #self.set_state(**current_state)

        likelihood_val = self._eval_raw()*self.norm 
        if likelihood_val < 0:
            print "Error: Pdf is 0 at state: ", current_state
            raise Exception("PdfVal")
        return likelihood_val


    def eval_data(self, data, **kwargs):
        """ Evaluate the pdf on a given value of data

        This is useful if one needs a function of data
        but doesn't know the name of the data parameter.
        Any other parameters can be set via kwargs, and
        the likelihood will be evaluated after those
        parameters are set.

        In the event that the kwargs set the data, the
        first argument of this function will over ride
        that setting.

        """
        self.set_state(**kwargs)
        self.set_data(data)
        return self.eval()


    def normalize(self):
        """ Integrate over the data
        at the current parameter point
        """

        # Check if the normalization has been cached
        # If so, return that cache
        norm_key = str(self.param_state().items())
        try:
            norm = self.normalization_cache[norm_key]
            self.norm = norm
            #self.logging.debug("Got Norm From Cache: %s from state: %s" % \
            #                       (self.norm, str(self.param_state())) )
            return norm
        except KeyError:
            pass

        # First we save the current value of data
        # since we are about to integrate over it
        current_data = self.get_data()
        
        # To normalize, we integrate the 'raw' function
        # over data
        def func_for_norm(data_val):
            self.set_data(data_val)
            return self._eval_raw()

        # If not, integrate over the data, invert it, 
        # and store the cache
        data_min, data_max = (self.data.min, self.data.max)
        #self.logging.debug("Integrating: " + str(self.param_state()))
        #integral, err = integrate.quad(self._eval_raw, data_min, data_max) 
        integral, err = integrate.quad(func_for_norm, data_min, data_max) 
        self.norm = 1.0 / integral

        if self.norm <= 0:
            print "Error: Normalization is <= 0"
            raise Exception("BadNormalization")

        #print "Found integral: ", integral, " Normalization: ", self.norm
        self.normalization_cache[norm_key] = self.norm

        # Restore the data value
        self.set_data(current_data)

        self.logging.debug("Got Norm From Integral: %s from state: %s" % \
                               (self.norm, str(self.param_state())) )

        return self.norm
    

    # Make the class callable:
    def __call__(self, *args, **kwargs):
        return self.eval(*args, **kwargs)


    # Log Likelihood
    def loglikelihood(self, *args, **kwargs):
        likelihood = self.eval(*args, **kwargs)
        if likelihood == 0:
            return -1*np.inf
            #return 0
        return log(likelihood)


    # Negative Log Likelihood
    def nll(self, *args, **kwargs):
        return -1*self.loglikelihood(*args, **kwargs)

    
    def get_interval(self, percentage):
        """ Get an inverval over the data parameter 
        which contains 'percentage' of the probability
        
        """
        points = self.data.linspace()
        pair_list = zip(points, map(self.eval_data, points))

        # ordering rule is maximum likelihood
        # sort by descending in likelihood
        pair_list = sorted(pair_list, key=lambda pair: pair[1], reverse=True)

        delta = (self.data.max - self.data.min) / self.data.num_points
        
        total_likelihood = 0.0
        accepted_point_list = []
        for (point, likelihood) in pair_list:
            accepted_point_list.append(point)
            total_likelihood += likelihood*delta
            if total_likelihood >= percentage: break
        
        interval = (min(accepted_point_list), max(accepted_point_list))
        return interval


    def get_neyman(self, percentage, param):
        """ Create a list of intervals
        for the parameter 'param'

        """

        param_var = self.param_dict[param] #getattr(self, param)
        
        interval_list = []
        for param_point in param_var.linspace():
            setattr(self, param, param_point)
            interval = self.get_interval(percentage)
            interval_list.append( (param_point, interval) )

        return interval_list


    def make_plot(self, interval=None):
        """ Plot the likelihood over data
        """

        #def eval_data(data):
        #    self.set_data(data)
        #    return self.eval()

        x = self.data.linspace()
        y = map(self.eval_data, x)
        plt.plot(x, y)
        plt.xlabel('x')
        plt.ylabel('likelihood(x)')
        x1,x2,y1,y2 = plt.axis()

        if interval != None:
            interval = self.get_interval(interval)
            plt.vlines(interval[0], y1, y2)
            plt.vlines(interval[1], y1, y2)

        return


    def invert_neyman(self, data, neyman=None, percentage=None, param=None, **kwargs):
        """ Invert neyman to get confidence interval
        
        The Neyman list looks like: [ (mu, (d0, d1)), ...
        """
        
        mu_list = []
        
        if neyman==None:
            neyman = self.get_neyman(percentage, param)

        for item in neyman:
            (mu, (d0, d1)) = item
            if d0 <= data and data <= d1:
                mu_list.append(mu)
            pass
        
        if len(mu_list)==0:
            print "Error: No acceptable values of param found"

        return (min(mu_list), max(mu_list))


    def fitTo(self, data, params, **kwargs):
        """ Minmize the supplied parameters based on the nll

        Set the values of the minimized parameters in the
        likelihood's 'state'.  Use any keyword arguments as
        initial values to parameters, or for any other
        (non-minimized) parameters in the model

        TO DO: Split the kwargs into args for the nll
        and args for optimize.minimize, and throw
        exceptions for all others
        (Should also warn about argument clashes with
        optimize when the likelihood function is initialized...)

        """

        # Log
        self.logging.debug( "Minimizing: " + str(params))

        # NOT YET IMPLEMENTED
        # Create a key based on the values of the params to not minimize
        # and the list of params to minimize (possibly overkill, but whatevs)
        # as well as the data being minimized
        constant_params = [item for item in self.param_state().items() if item[0] not in params]
        #cache_key = hash( (data, frozenset(constant_params), frozenset(params)) )
        cache_key = (data, frozenset(constant_params), frozenset(params))
        if cache_key in self.minimization_cache:
            state = self.minimization_cache[cache_key] 
            self.set_state(**state)
            return

        # Minimize the supplied params
        if len(params)==0:
            print "Error: No Paramaterize to Minimize"
            raise Exception("FitTo")
            return

        # Set the value of the data to the supplied data
        self.set_data(data)

        current_state = self.total_state()
        self.normalize()

        # Create the function for minimization
        def nnl_for_min(param_value_list):
            """ Create the wrapper function for scipy.optimize

            """

            for (param, val) in zip(params, param_value_list):
                setattr(self, param, val)

            # nll without normalization
            return -1*log(self._eval_raw())

        # Get the initial guess
        guess = [getattr(self, param) for param in params]

        # Run the minimization
        res = optimize.minimize(nnl_for_min, guess)
        self.logging.debug("Successfully Minimized the function to value:", res)
        
        # Set the values to the minimum
        min_values = res.x
        for (param, val) in zip(params, min_values):
            self.logging.debug("Minimized value of %s : %s" % (param, val))
            setattr(self, param, val)

        # Cache the result
        self.minimization_cache[cache_key] = self.total_state()

        return


    def profile(self, poi, nuisance, **kwargs):
        """ Return the profile likelihood as a function of the poi

        (parameter of interest), minimizing over the nuisance parameters

        return the nll of the profile likelihood
        """

        self.set_state(**kwargs)

        # Save the current value since we are evaluating
        # the profile likelihood at this point
        current_poi_value = getattr(self, poi)

        # Get the set of parameters
        all_params = [poi]
        all_params.extend(nuisance)
        self.logging.debug( "Profiling poi: %s and nuisance params %s" % (poi, str(nuisance)) )

        # Save the current state
        saved_state = {}
        for arg in all_params:
            saved_state[arg] = getattr(self, arg)

        # Get the constant parameters
        '''
        const_params = []
        for param in self._arg_list:
            if param == poi: continue
            if param in nuisance: continue
            const_params.append( (param, getattr(self, param)) )
        const_params = tuple(const_params)
        '''

        # Get the global min
        #cache_key = hash((self.get_data(), frozenset(all_params))) 
        cache_key = (self.get_data(), frozenset(all_params)) 
        if cache_key in self.nll_cache:
            global_nll = self.nll_cache[cache_key]
        else:
            self.fitTo(self.get_data(), params=all_params)
            global_nll = self.nll()
            self.nll_cache[cache_key] = global_nll
        
        # Get the local min at this point
        setattr(self, poi, current_poi_value)
        self.fitTo(self.get_data(), params=nuisance)
        local_nll = self.nll()

        # Restore the state
        for arg in all_params:
            setattr(self, arg, saved_state[arg])

        '''
        nll = self.nll(data)
        output = "Profile: global nll: " + str(global_nll) \
            + " local_nll: " + str(local_nll) \
            + " original nll: " + str(nll)
        self.logging.debug(output)
        '''
        self.logging.debug("Found Profile. Local nll: %s Global nll: %s" % (local_nll, global_nll))

        return local_nll - global_nll


    def sample(self, params, nsamples=1, method='mc', **kwargs):
        """ Generate sample points based on the likelihood

        """

        self.set_state(**kwargs)

        if method=='':
            print "Must supply method for generate()"
            raise Exception("generate")
        elif method=='mcmc':
            return self.sample_mcmc(params, nsamples, **kwargs)
        elif method=='mc':
            return self.sample_mc(params, nsamples, **kwargs)
        else:
            print "Supplied invalid method for generate(): %s" % method
            raise Exception("generate")


    def sample_mc(self, params=[], nsamples=1, **kwargs):
        """ Sample from the likelihood using brute force Monte-Carlo

        """

        self.set_state(**kwargs)

        for param in params:
            if not hasattr(self, param):
                self.log.error("Cannot sample parameter %s, no such parameter" % param)
                raise Exception("SampleError")
            pass

        # Save the current state
        saved_state = self.total_state()

        results=[]
        for i_sample in xrange(nsamples):

            # Some output
            if i_sample % 1000 == 1:
                print "Generating sample: %s" % i_sample

            while True:
            
                # Set the values
                for param in params:
                    param_var = self.args[param] #param_dict[param]
                    #(param_min, param_max) = (param_var.
                    val = random.uniform(param_var.min, param_var.max)
                    #self.logging.debug("Setting attribute: %s %s" % (param, val))
                    setattr(self, param, val)
            
                # Get the likelihood
                lhood = self.eval()

                # Throw the Monte-Carlo dice:
                mc_val = random.uniform(0.0, 1.0)

                '''
                debug_string =  "MC Accept/Reject: "
                total_state = self.total_state()
                debug_string += " state: " + str(total_state)
                debug_string += " likelihood: %s" % lhood
                debug_string += " mc_val: %s" % mc_val
                self.logging.debug(debug_string)
                '''

                if lhood > mc_val: break

            point = {}
            for param in params:
                point[param] = getattr(self, param)

            results.append(point)

        self.set_state(**saved_state)
        return results


    def sample_mcmc(self, params=[], nsamples=1, nwalkers=6, burn_in=500):
        """
        
        Generate 'nsamples' points based on the likelihood
        Use the emcee package for MarkovChain Monte-Carlo

        return a list of dictionaries of name, val pairs for the
        supplied points
        """

        try:
            import emcee
        except ImportError:
            print "Cannot use Markov Chain Monte Carlo, must install 'emcee' package"
            print "See: https://danfm.ca/emcee/"
            raise

        saved_state = self.total_state()

        def func_for_emcee(val_list):
            """ Requires the log probability and 
            params as an array """
            
            # Set the state based on the input list
            for (param, val) in zip(params, val_list):
                setattr(self, param, val)
                
            # Return the log likelihood
            return self.loglikelihood()

        # Setup emcee
        # WARNING: Check ranges here
        # Set up the initial states of the walkers:
        p0 = [[random.uniform(param.min, param.max) for (name, param) in self.args.iteritems() if name in params] for j in xrange(nwalkers)]
        print "Initial guess: ", p0
        print self.param_dict
        #p0 = [[random.uniform(-1, 1) for i in params] for j in xrange(nwalkers)]
    
        #nwalkers
        ndim = len(params)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, func_for_emcee)

        # Run 500 steps as a burn-in.
        pos, prob, state = sampler.run_mcmc(p0, burn_in)

        # Reset the chain to remove the burn-in samples.
        sampler.reset()
        
        # Starting from the final position in the burn-in chain, sample for 2000
        # steps.
        sampler.run_mcmc(pos, nsamples, rstate0=state)
        
        # Unpack the results
        results = []
        for (isample,sample) in enumerate(sampler.flatchain):
            point = {}
            for (iparam, param) in enumerate(params):
                point[param] = sample[iparam]
            results.append(point)
        
        self.set_state(**saved_state)
        return results
        

