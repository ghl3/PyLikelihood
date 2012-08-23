
import sys
import math

import logging
import inspect
from pprint import pprint

import scipy.optimize
import scipy.integrate

import random


class Likelihood(object):

    _arg_list=[]
    _likelihood_function=None
    _integral_value=1.0
    _var_ranges = {}
    _cache = {}
    log = logging.getLogger()

    def __init__(self, func=None):
        if func == None:
            return
        self.setLikelihood(func)

    def setRange(self, param, min, max):
        if param in self._arg_list:
            self._var_ranges[param] = (min, max)
        else:
            self.log.error("Cannot set range, param %s not found" % param)
            raise Exception("ParamNotFound")

    #def ActivateLogging(level=logging.DEBUG):
    #    FORMAT = "%(message)s"
    #    logging.basicConfig(format=FORMAT)
    #    logging.root.setLevel(logging.DEBUG)

    def setLikelihood(self, func, **kwargs):
        """ Set the likelihood function to be 'func'

        The first argument is interpreted as a list
        that represents the data
        
        This instance of the class will get members for
        every parameter of the function (except for the
        first one).
        Future calls to the method "likelihood" will call
        this function, but using the stored values in this
        class as the arguments.

        keyword arguments set the default values to
        these parameters.  Otherwise, they get None

        """
        
        # Get the parameters of the function
        func_spec = inspect.getargspec(func)
        arg_list = func_spec.args[1:]
        for arg_name in arg_list:
            if arg_name in kwargs:
                setattr(self, arg_name, kwargs[arg_name])
            else:
                setattr(self, arg_name, None)

        for arg_name in kwargs:
            if arg_name not in arg_list:
                self._log.error("Error: setLikelihood recieved argument %s" % arg_name /
                               " but this is not a keyword argument of the likelihood" )
                raise Exception("Likelihood Argument Error")

        self._arg_list = arg_list
        self._likelihood_function = func


    def get_function_args(self):
        """ Return the kwargs for the callable
        based on the current state
        """
        function_args = {}
        for arg in self._arg_list:
            val = getattr(self, arg)
            function_args[arg] = val
        return function_args


    def set_state(self, **kwargs):
        """ Set the state based on the values of the arguments

        Return any args that aren't parameters of the likelihood
        """
        unused_args = {}
        for (arg, val) in kwargs.iteritems():
            if arg in self._arg_list:
                setattr(self, arg, val)
            else:
                unused_args[arg] = val
        return unused_args


    def get_state(self):
        state = {}
        for arg in self._arg_list:
            state[arg] = getattr(self, arg)
        return state        

    def eval(self, data_point, **kwargs):
        """ Evaluate the callable function on a single dataset
        
        """

        # Set the state based on the kwargs
        self.set_state(**kwargs)

        # Get the value
        func_args = self.get_function_args()
        likelihood_val = self._likelihood_function(data_point, **func_args)

        self.check_value(likelihood_val)

        return likelihood_val


    def likelihood(self, dataset, **kwargs):
        """ Call the likelihood function on the supplied dataset

        Dataset should be a list of data points

        Use the callable function that was previously set, 
        and use the current 'state' of the class for the
        arguments to that function.

        Any supplied keyword arguments that match arguments
        to that function will be set as arguments to the 
        function, and they will be saved in the state

        """

        self.set_state(**kwargs)

        # Get the value
        func_args = self.get_function_args()
        likelihood_val = 1.0
        for point in dataset:
            likelihood_val *= self._likelihood_function(point, **func_args)
            
        self.check_value(likelihood_val)
        
        return likelihood_val


    def loglikelihood(self, dataset, **kwargs):
        return math.log(self.likelihood(dataset, **kwargs))

    def nll(self, dataset, **kwargs):
        """ Return the negative log likelihood 

        Evaluate on the supplied data and return the value
        We catch any ValueError exceptions here since these
        can occur 

        """

        self.set_state(**kwargs)

        likelihood_val = self.likelihood(dataset, **kwargs)

        try:
            neg_log_val = -1*math.log(likelihood_val)
        except ValueError:
            self.log.error("Encountered Val Error.  Input to log: ", likelihood_val)
            raise Exception("NegativeLogLikelihoodEval")

        return neg_log_val


    def minimize(self, dataset, params, **kwargs):
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

        # Set the current state based on keyword args
        # This will effect the initial guess and any
        # constant parameters (non-minimized)
        unused_params = self.set_state(**kwargs)

        if len(params)==0:
            return

        # Create the function for minimization
        def nnl_for_min(param_values):
            """ Create the wrapper function for scipy.optimize

            """

            for (nuis, val) in zip(params, param_values):
                setattr(self, nuis, val)
            #print self.get_state()
            return self.nll(dataset)

        # Get the initial guess
        guess = [getattr(self, param) for param in params]
        self.log.debug("Minimizing: ")
        for param, val in zip(params, guess):
            self.log.debug("%s : %s" % (param, val))

        # Run the minimization
        res = scipy.optimize.minimize(nnl_for_min, guess, **unused_params)
        self.log.debug("Successfully Minimized:", res)

        #bounds = [self.bounds[param] for param in params]
        #print "Minimizing: ", zip(params, guess, bounds)
        #res = scipy.optimize.minimize(nnl_for_min, guess, bounds=bounds, method='SLSQP')
        
        # Set the values to the minimum
        min_values = res.x
        for (param, val) in zip(params, min_values):
            self.log.debug("%s : %s" % (param, val))
            setattr(self, param, val)

        return min


    def profile(self, dataset, poi, nuisance, **kwargs):
        """ Return the profile likelihood as a function of the poi
        (parameter of interest), minimizing over the nuisance parameters

        return the nll of the profile likelihood
        """

        #if len(nuisance)==0:
        #    print "Error: Must supply nuisance parameters"
        #    raise Exception("ProfileLikelihood")

        # Set the State
        self.set_state(**kwargs)

        # Save the current value
        current_poi_value = getattr(self, poi)

        # Get the set of parameters
        all_params = [poi]
        all_params.extend(nuisance)
        self.log.debug( "All Params: %s" % all_params)

        # Get the constant parameters
        const_params = []
        for (param) in self._arg_list:
            if param == poi: continue
            if param in nuisance: continue
            const_params.append( (param, getattr(self, param)) )
        const_params = tuple(const_params)

        # Get the global min
        if const_params not in self._cache:
            global_min = self.minimize(dataset, params=all_params, **kwargs)
            global_nll = self.nll(dataset, **kwargs)
            self._cache[const_params] = global_nll
        else:
            global_nll = self._cache[const_params]

        # Get the local min at this point
        setattr(self, poi, current_poi_value)
        local_min = self.minimize(dataset, params=nuisance, **kwargs)
        local_nll = self.nll(dataset, **kwargs)

        return local_nll #- global_nll


    def integral(self, dataset, param, range=None):
        """ Get the integral of the function evaluated over data

        """
        
        if range == None and param not in self._var_ranges:
            self.log.error("Integral: Param %s has no range" % param)
            raise Exception("IntegralRange")

        if range==None:
            range = self._var_ranges[param]

        saved_state = self.get_state()

        def func_for_integral(val):
            setattr(self, param, val)
            return self.eval(dataset)

        result = scipy.integrate.quad(func_for_integral, range[0], range[1])

        self.set_state(**saved_state)

        return result[0]


    def sample(self, dataset, params, nsamples=1, method='mc', **kwargs):
        """
        Generate sample points based on the given likelihood

        """

        if method=='':
            print "Must supply method for generate()"
            raise Exception("generate")
        elif method=='mcmc':
            return self.sample_mcmc(dataset, params, nsamples, **kwargs)
        elif method=='mc':
            return self.sample_mc(dataset, params, nsamples, **kwargs)
        else:
            print "Supplied invalid method for generate(): %s" % method
            raise Exception("generate")

    def sample_mcmc(self, dataset, params=[], nsamples=1, nwalkers=6):
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

        saved_state = self.get_state()

        def func_for_emcee(val_list):
            """ Requires the log probability and 
            params as an array """
            
            # Set the state based on the input list
            for (param, val) in zip(params, val_list):
                setattr(self, param, val)
                
            # Return the log likelihood
            log_lhood = self.loglikelihood(dataset)
            return log_lhood

        # Setup emcee
        # WARNING: Check ranges here
        p0 = [[random.uniform(-1, 1) for i in params] for j in xrange(nwalkers)]
    
        #nwalkers
        ndim = len(params)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, func_for_emcee)

        # Run 500 steps as a burn-in.
        pos, prob, state = sampler.run_mcmc(p0, 500)

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
        

    def sample_mc(self, dataset, params=[], nsamples=1):
        
        for param in params:
            if param not in self._var_ranges:
                self.log.error("Cannot sample parameter %s, must supply range" % param)
                raise Exception("SampleError")
            pass

        # Save the current state
        saved_state = self.get_state()

        results=[]
        for i_sample in xrange(nsamples):
            while True:
            
                # Set the values
                for param in params:
                    range = self._var_ranges[param]
                    val = random.uniform(range[0], range[1])
                    setattr(self, param, val)
            
                # Get the likelihood
                lhood = self.likelihood(dataset)

                # Throw the Monte-Carlo dice:
                mc_val = random.uniform(0.0, 1.0)

                if lhood > mc_val:
                    break

            point = {}
            for param in params:
                point[param] = getattr(self, param)

            results.append(point)

        self.set_state(**saved_state)
        return results


    def check_value(self, val):
        if val <= 0.0:
            self.log.error("Error: Likelihood evaluated to < 0: ", val)
            raise Exception("LikelihoodEval")
        if math.isnan(val):
            self.log.error("Error: Likelihood value is NAN: ", val)
            raise Exception("LikelihoodEval")
        if math.isinf(val):
            self.log.error("Error: Likelihood value is INF: ", val)
            raise Exception("LikelihoodEval")
        return

