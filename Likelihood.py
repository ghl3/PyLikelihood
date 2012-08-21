
import sys
import math

import inspect
from pprint import pprint

import scipy.optimize
#from scipy import stats
from scipy.stats import poisson
from scipy.stats import norm

import numpy as np

def pois(x, l):
    rv = poisson([l])
    return rv.pmf(x)[0]

def gauss(x, mu, sigma):
    rv = norm(loc=mu, scale=sigma)
    return rv.pdf([x])[0]

class Likelihood(object):
    def __init__(self):
        self._arg_list=[]
        self._likelihood_function=None
        #self.n=None
        #self.mu=None
        #self.delta=None
        #self.alpha=None

        #self.bounds={}
        #self.bounds["n"] = (0, None)
        #self.bounds["mu"] = (0, 10)
        #self.bounds["alpha"] = (-5, 5)
        #self.bounds["delta"] = (None, None)
        
    def print_state(self):
        pprint(vars(self))

        
    def nll(self, dataset, **kwargs):
        val = 0.0
        for point in dataset:
            try:
                val += -1*math.log(self.likelihood(point, **kwargs))
            except ValueError:
                return 999999 #np.inf # Inf #sys.float_info.max
        return val


    def minimize(self, dataset, params=[], **kwargs):

        # Create the function for minimization
        def nnl_for_min(param_values):
            """ Minimize the likelihood over all supplied params

            """
            
            # Set the value of the var and the nuisance
            for (nuis, val) in zip(params, param_values):
                setattr(self, nuis, val)
            self.print_state()
            return self.nll(dataset, **kwargs)

        # Get the global minimum
        guess = [getattr(self, param) for param in params]
        print "Minimizing: ", zip(params, guess)
        res = scipy.optimize.minimize(nnl_for_min, guess)

        #bounds = [self.bounds[param] for param in params]
        #print "Minimizing: ", zip(params, guess, bounds)
        #res = scipy.optimize.minimize(nnl_for_min, guess, bounds=bounds, method='SLSQP')
        
        # Set the values to the minimum
        min_values = res.x
        
        print "Successfully Minimized:", res
        for (param, val) in zip(params, min_values):
            print param, val
            setattr(self, param, val)

        return min

    #
    # Experimental Methods
    #


    def SetLikelihood(self, func):
        """ Set the likelihood function
        
        This instance of the class will get members for
        every parameter of the function, and future calls
        to the likelihood will use the current values
        stored in this instance

        The first argument is interpreted as a list
        that represents the data

        """
        
        # Get the parameters of the function
        func_spec = inspect.getargspec(func)
        arg_list = func_spec.args[1:]
        for arg in arg_list:
            setattr(self, arg, None)

        self._arg_list = arg_list
        self._likelihood_function = func


    def likelihood(self, data, **kwargs):
        """ Call the likelihood function with the
        current state of the class

        """

        # Take any keyword arguments to CallLikelihood
        # and set the current value of the class before
        # actually calling the likelihood
        for (arg, val) in kwargs.iteritems():
            if arg in self._arg_list:
                setattr(self, arg, val)
            else:
                print "Error: %s is not a parameter of the likelihood" % arg
                raise Exception("Bad Likelihood Parameter")

        # Create the argumets
        kw_args = {}
        for arg in self._arg_list:
            val = getattr(self, arg)
            kw_args[arg] = val
        return self._likelihood_function(data, **kw_args)


        
'''
    def profile(self, dataset, var, nuisance=[]):
        """ Return the profile likelihood of 'var',
        minimizing over all parameters in 'profile'
        """
'''
