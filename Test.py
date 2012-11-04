#!/usr/bin/env python

import math
import itertools
import operator

from Likelihood import *

import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import poisson
from scipy.stats import norm
import scipy.integrate

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
    small_num =  10e-25
    if math.isnan(val):
        #print "Val is NAN"
        return small_num
    if val <= small_num:
        #print "Val is small", val
        return small_num

    return val

def my_func(x, y):
    return math.exp(-1*(x-y)*(x-y))
    #return gauss(x, mu, sigma)

def triple_gauss(x, y, z, mu=0, sigma=1):
    return gauss(x,mu,sigma)*gauss(y,mu,sigma)*gauss(z,mu,sigma)


def makeGrid(num_points, ranges):
    """ Create a list of dimension 'dim' points on the unit cube
    
    dim = dimension of each vector (ie, the number of alpha's)
    num_points = the number of points in each direction
    min, max are the range of the box
    
    ranges  [(min1, max1), (min2, max2)...]

    Each dimension has num_points that are equally separated
    (0, 0, .2), (0, 0, .4), ...
    (0, .2, 0), (0, .4, 0), ...
    ...
    
    ps. itertools rocks
    
    """

    # Make the list of deltas
    deltas = [ (pair[1]-pair[0])/float(num_points-1) for pair in ranges]

    # Create the iterator factory
    def createGenerator(dx):
        for j in xrange(num_points):
            yield j*dx

    # Make the list of iteratoers
    iteraters = [createGenerator(dx) for dx in deltas]

    # return itertools magic
    return itertools.product( *iteraters )


def main():


    params = ["mu", "alpha_1", "alpha_2"]
    npoints = 10
    ranges = [(-5,5), (-5,5), (-5,5)]
    volume = reduce(operator.mul, [(pair[1]-pair[0])/float(npoints) for pair in ranges] )
    grid = makeGrid(npoints, ranges)

    integral = 0.0
    for point in grid:
        #val = simple_likelihood(10, 8, point[0], point[1], 1, point[2], 2) 
        val = triple_gauss(*point)
        integral += val*volume

    print integral
    return
                   

    # Try normalizing the 'simple_likelihood'
    #params = ["mu", "alpha_1", "alpha_2"]
    #grid = makeGrid( len(params), .1, 

    

    FORMAT = "%(message)s"
    logging.basicConfig(format=FORMAT)
    logging.root.setLevel(logging.DEBUG)

    # Create the dataset
    dataset = [2]

    # Test the integral
    # Create a simple likelihood model

    model = Likelihood(dataset, gauss)
    model.mu=0.0
    model.sigma=1.0
    print model.get_state()

    print model.eval(2)
    print model.eval(2, mu=3)

    print model.likelihood()
    model.minimize(params=['mu'])

    # Sample the model:
    model.setRange('mu', -5, 5)

    print "Scipy Integral: ", model.integral(2, 'mu')

    # Test the sampling
    nsamples=500
    samples = model.sample(params=['mu'], nsamples=nsamples, method='mcmc')
    lik_samples = [model.likelihood(**sample) for sample in samples]
    mu_samples = [sample['mu'] for sample in samples]

    #for i in range(5000):
    #    sample = model.sample(dataset, args=['mu'])
    #    print sample
    #    likelihood_samples.append(sample['likelihood'])
    #    mu_samples.append(sample['values']['mu'])

    plt.figure()
    # Get the "true" values
    points = scipy.linspace(-5, 5, num=100)
    y = [model.likelihood(mu=p) for p in points]
    plt.plot(points, y, label="likelihood")
    # Plot the histogram of samples
    plt.hist(mu_samples, 100, normed=True, facecolor='g', alpha=0.75, label="sampled")
    plt.legend()
    plt.savefig("sample_test.pdf")    

    return

    #
    # Create a more complicated likelihood model
    #

    model = Likelihood(simple_likelihood)

    model.n = 10
    model.mu = 1.0
    model.alpha_1 = 0
    model.delta_1 = .2
    model.alpha_2 = 0
    model.delta_2 = .3

    data = [12]

    # Test the integral
    '''
    integral = 0.0
    for point in range(0, 100):
        integral += model.eval(point)
    print "Integral: ", integral
    print "Scipy Integral: ", model.integral()
    return
    '''



    pll = model.profile(data, "mu", nuisance=['alpha_1', 'alpha_2'])
    print "Profile Likelihood: ", pll

    print model.likelihood(data)
    print model.get_state()

    print model.likelihood(data, alpha_1=1)
    print model.get_state()

    model.alpha_1=0.0

    # Test the minimization
    #model.minimize(data, params=['mu','alpha'])
    model.minimize(data, params=['mu'])
    model.minimize(data, params=['mu','alpha_1'])
    model.minimize(data, params=['mu','alpha_1', 'alpha_2'])
                    
    # Plot the likelihood as a function of mu
    x = scipy.linspace(0, 2, num=100)
    #y = [model.profile(data, mu=p) for p in x]
    model.log.setLevel(logging.WARNING)
    y = [model.nll(data, mu=p) for p in x]
    z = [model.profile(data, "mu", ["alpha_1", "alpha_2"], mu=p) for p in x]
    model.log.setLevel(logging.DEBUG)
    #z = [model.profile(data, "mu", ['alpha_1'], mu=p) for p in x]
    plt.figure()
    plt.plot(x,y, label="nll")
    plt.plot(x,z, label="profile")
    plt.legend()
    plt.savefig("nll.pdf")    
    return

    # Test the pois
    x = np.arange(-5, 5, .1)
    y = [gauss(point, 0, 1) for point in x]
    plt.figure()
    plt.plot(x, y)
    plt.savefig("gauss.pdf")


    #x = scipy.linspace(0,10,11)
    #pylab.plot(x, pois(x,5))

if __name__=="__main__":
    main()
