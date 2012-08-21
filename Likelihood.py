
import sys
import math

import inspect
from pprint import pprint

import scipy.optimize

class Likelihood(object):
    def __init__(self, func=None):
        self._arg_list=[]
        self._likelihood_function=None
        if func == None:
            return
        self.SetLikelihood(func)

    def print_state(self):
        pprint(vars(self))

    def SetLikelihood(self, func, **kwargs):
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
                print "Error: SetLikelihood recieved argument %s" % arg_name,
                print " but this is not a keyword argument of the likelihood"
                raise Exception("Likelihood Argument Error")

        self._arg_list = arg_list
        self._likelihood_function = func


    def likelihood(self, data, **kwargs):
        """ Call the likelihood function on the supplied data

        Use the callable function that was previously set, 
        and use the current 'state' of the class for the
        arguments to that function.

        Any supplied keyword arguments that match arguments
        to that function will be set as arguments to the 
        function, and they will be saved in the state

        """

        # Set the state based on keyword arguments
        for (arg, val) in kwargs.iteritems():
            if arg in self._arg_list:
                setattr(self, arg, val)
            else:
                print "Error: %s is not a parameter of the likelihood" % arg
                raise Exception("Bad Likelihood Parameter")

        # Create the argumets to the function
        kw_args = {}
        for arg in self._arg_list:
            val = getattr(self, arg)
            kw_args[arg] = val

        likelihood_val = self._likelihood_function(data, **kw_args)

        if likelihood_val < 0:
            print "Error: Likelihood evaluated to < 0: ", likelihood_val
            raise Exception("LikelihoodEval")
        if math.isnan(likelihood_val):
            print "Error: Likelihood value is NAN: ", likelihood_val
            raise Exception("LikelihoodEval")
        if math.isinf(likelihood_val):
            print "Error: Likelihood value is INF: ", likelihood_val
            raise Exception("LikelihoodEval")

        return likelihood_val


    def nll(self, data, **kwargs):
        """ Return the negative log likelihood 

        Evaluate on the supplied data and return the value
        We catch any ValueError exceptions here since these
        can occur 

        """

        val = 0.0
        for point in data:
            likelihood_val = self.likelihood(point, **kwargs)
            if likelihood_val == 0.0:
                print "Error: Likelihood evaluates to 0.0"
                raise Exception("LikelihoodEval")
            try:
                val += -1*math.log(likelihood_val)
            except ValueError:
                print "Encountered Val Error.  Input to log: ", likelihood_val
                raise Exception("NegativeLogLikelihoodEval")
        return val


    def minimize(self, data, params, **kwargs):
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

        print "Minimizing the following parameters:"
        params

        # Create the function for minimization
        def nnl_for_min(param_values):
            """ Create the wrapper function for scipy.optimize

            """

            for (nuis, val) in zip(params, param_values):
                setattr(self, nuis, val)
            self.print_state()
            return self.nll(data, **kwargs)

        # Get the initial guess
        guess = [getattr(self, param) for param in params]
        print "Minimizing: ", zip(params, guess)

        # Run the minimization
        res = scipy.optimize.minimize(nnl_for_min, guess)
        print "Successfully Minimized:", res
        #bounds = [self.bounds[param] for param in params]
        #print "Minimizing: ", zip(params, guess, bounds)
        #res = scipy.optimize.minimize(nnl_for_min, guess, bounds=bounds, method='SLSQP')
        
        # Set the values to the minimum
        min_values = res.x
        for (param, val) in zip(params, min_values):
            print param, val
            setattr(self, param, val)

        return min


    def profile(self, data, poi, nuisance, **kwargs):
        """ Return the profile likelihood as a function of the poi
        (parameter of interest), minimizing over the nuisance parameters

        return the nll of the profile likelihood
        """

        if len(nuisance)==0:
            print "Error: Must supply nuisance parameters"
            raise Exception("ProfileLikelihood")

        current_poi_value = getattr(self, poi)

        # Get the global minimum
        all_params = [poi]
        all_params.extend(nuisance)
        print "All Params: ", all_params
        global_min = self.minimize(data, params=all_params, **kwargs)
        global_nll = self.nll(data)

        # Get the local min at this point
        setattr(self, poi, current_poi_value)
        local_min = self.minimize(data, params=nuisance, **kwargs)
        local_nll = self.nll(data)

        return local_nll - global_nll
