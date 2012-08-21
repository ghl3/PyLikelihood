#!/usr/bin/env python

import math

from Likelihood import *

import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import poisson
from scipy.stats import norm


def pois(x, l):
    rv = poisson([l])
    return rv.pmf(x)[0]

def gauss(x, mu, sigma):
    rv = norm(loc=mu, scale=sigma)
    return rv.pdf([x])[0]

def simple_likelihood(d, n, mu, alpha_1, delta_1, alpha_2, delta_2):
    """ The probability of a single data point given parametres
    
    """
    n_hat = n*mu*(1.0 + alpha_1*delta_1 + alpha_2*delta_2)
    val = pois(d, n_hat)*gauss(0.0, alpha_1, 1.0)*gauss(0.0, alpha_2, 1.0)

    # Add some simple protection
    small_num = .00000000000000000000000000000000001 
    if math.isnan(val):
        print "Val is NAN"
        return small_num
    if val <= small_num:
        print "Val is small", val
        return small_num

    return val


def my_func(x, mu, sigma):
    return gauss(x, mu, sigma)


def main():

    model = Likelihood(my_func)
    model.mu=0.0
    model.sigma=1.0
    model.print_state()
    #model.bracket['sigma'] = (0, None)    
    dataset = [2]

    print model.likelihood(dataset)
    model.minimize(dataset, params=['mu'])

    # Test the pois
    '''
    x = np.arange(-5, 5, .1)
    y = [gauss(point, 0, 1) for point in x]
    plt.figure()
    plt.plot(x, y)
    plt.savefig("gauss.pdf")
    '''

    # Create a likelihood model
    #model = Likelihood()
    #model.SetLikelihood(simple_likelihood)
    model = Likelihood(simple_likelihood)

    model.n = 10
    model.mu = 1.0
    model.alpha_1 = 0
    model.delta_1 = 2
    model.alpha_2 = 0
    model.delta_2 = 3

    data = [12]

    pll = model.profile(data, "mu", nuisance=['alpha_1', 'alpha_2'])
    print "Profile Likelihood: ", pll

    print model.likelihood(data)
    model.print_state()

    print model.likelihood(data, alpha_1=1)
    model.print_state()

    model.alpha_1=0.0

    # Test the minimization
    #model.minimize(data, params=['mu','alpha'])
    model.minimize(data, params=['mu'])
    model.minimize(data, params=['mu','alpha_1'])
    model.minimize(data, params=['mu','alpha_1', 'alpha_2'])
    return
                    
    # Plot the likelihood as a function of mu
    x = scipy.linspace(0, 2, num=100)
    y = [model.nll(data, mu=p) for p in x]
    #pylab.plot(x, y)
    plt.figure()
    plt.plot(x,y)
    #plt.axis([0, 3, 0, 20])
    plt.savefig("nll.pdf")    
    return


    #x = scipy.linspace(0,10,11)
    #pylab.plot(x, pois(x,5))

if __name__=="__main__":
    main()
