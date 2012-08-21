#!/usr/bin/env python

import math

import scipy
from scipy import stats
from scipy.stats import poisson
from scipy.stats import norm

import numpy as np
import matplotlib.pyplot as plt

import pylab


def pois(x, l):
    rv = poisson([l])
    return rv.pmf(x)

def gauss(x, mu, sigma):
    rv = norm(loc=mu, scale=sigma)
    return rv.pdf([x])

class nll(object):
    def __init__(self):
        self.n=None
        self.mu=None
        self.delta=None
        self.alpha=None

    def _likelihood(self, d, **kwargs):
        """ The probability of a single data point given parametres

        """

        # Set the values of any given args
        for (member, val) in kwargs.iteritems():
            setattr(self, member, val)
        
        n_hat = self.n*self.mu*(1.0 + self.alpha*self.delta)
        return pois(self.n, n_hat)*gauss(self.alpha, 0.0, 1.0)
        
    def eval(self, dataset, **kwargs):
        val = 0.0
        for point in dataset:
            val += -1* math.log(self._likelihood(point, **kwargs))
        return val


    def minimize(self, dataset, params=[]):

        # Create the function for minimization
        def nnl_for_min(param_values):
            """ Minimize the likelihood over all supplied params

            """
            
            # Set the value of the var and the nuisance
            for (nuis, val) in zip(params, param_values):
                setattr(self, nuis, val)

            return self.eval(dataset)

        # Get the global minimum
        guess = [getattr(self, param) for param in params]
        print "Guess: ", guess
        res = scipy.optimize.minimize(nnl_for_min, guess)
        
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



def main():

    my_nll = nll()
    my_nll.n = 100
    my_nll.mu = 1.0
    my_nll.alpha = 0
    my_nll.delta = 2

    data = [110]
   
    # Test the minimization
    my_nll.minimize(data, params=['alpha', 'mu'])

    return

                    
 
                    
    # Plot the likelihood as a function of mu
    x = scipy.linspace(0,2,num=100)
    y = [my_nll.eval(data, mu=p) for p in x]
    pylab.plot(x, y)
    plt.savefig("bob.pdf")    
    return


    x = np.arange(0, 10)
    #x = np.arange(0, np.minimum(rv.dist.b, 3))
    plt.plot(x, pois(x,4))
    plt.savefig("frank.pdf")

    #x = scipy.linspace(0,10,11)
    #pylab.plot(x, pois(x,5))

    
if __name__ == "__main__":
    main()
    
