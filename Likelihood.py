
import sys
import math

import scipy.optimize
#from scipy import stats
from scipy.stats import poisson
from scipy.stats import norm

import numpy as np

from pprint import pprint


def pois(x, l):
    rv = poisson([l])
    return rv.pmf(x)[0]

def gauss(x, mu, sigma):
    rv = norm(loc=mu, scale=sigma)
    return rv.pdf([x])[0]

class Likelihood(object):
    def __init__(self):
        self.n=None
        self.mu=None
        self.delta=None
        self.alpha=None

        #self.bounds={}
        #self.bounds["n"] = (0, None)
        #self.bounds["mu"] = (0, 10)
        #self.bounds["alpha"] = (-5, 5)
        #self.bounds["delta"] = (None, None)
        
    def print_state(self):
        pprint(vars(self))

    def _likelihood(self, d, **kwargs):
        """ The probability of a single data point given parametres

        """
        # Set the values of any given args
        for (member, val) in kwargs.iteritems():
            setattr(self, member, val)
        
        n_hat = self.n*self.mu*(1.0 + self.alpha*self.delta)
        return pois(d, n_hat)*gauss(0.0, self.alpha, 1.0)

        
    def nll(self, dataset, **kwargs):
        val = 0.0
        for point in dataset:
            try:
                val += -1*math.log(self._likelihood(point, **kwargs))
            except ValueError:
                return 999999 #np.inf # Inf #sys.float_info.max
        return val


    def minimize(self, dataset, params=[]):

        # Create the function for minimization
        def nnl_for_min(param_values):
            """ Minimize the likelihood over all supplied params

            """
            
            # Set the value of the var and the nuisance
            for (nuis, val) in zip(params, param_values):
                setattr(self, nuis, val)

            self.print_state()

            return self.nll(dataset)

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
        
        
'''
    def profile(self, dataset, var, nuisance=[]):
        """ Return the profile likelihood of 'var',
        minimizing over all parameters in 'profile'
        """
'''
