#!/usr/bin/env python

from __future__ import division

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import poisson
from scipy.stats import norm
from scipy import integrate

from math import log

class variable:
    def __init__(self, name, min, max, num_points=1000):
        self.name = name
        self.min = min
        self.max = max
        self.num_points = num_points

    def linspace(self):
        return np.linspace(self.min, self.max, self.num_points)
    pass


class likelihood:
    def __init__(self, func, data_var):
        self.func = func
        self.data_var = data_var

        # Get the list of parameters
        func_spec = inspect.getargspec(pdf)
        self.param_list = [arg for arg in func_spec.args if arg != data_var]
        
    def normalize():
    
   def eval(x):
       


        self.max = max
        self.num_points = num_points
    

# Annoying boilerplate
def pois(x, mu):
    return poisson.pmf(x, mu)

def gauss(x, mu, sigma):
    return norm.pdf(x, mu, sigma)




# Create the likelihood
def likelihood(d, mu=5.0, mu0=5.0, sigma=2.0):
    #return pois(d, mu)*gauss(mu0, mu, sigma)
    return likelihood.norm*gauss(d, mu, 2.0)*gauss(mu0, mu, sigma)
likelihood.norm = 1.0
#likelihood.cache = {}

likelihood_int, err = integrate.quad(likelihood, -np.inf, np.inf) 
likelihood.norm = 1.0 / likelihood_int
print likelihood.norm

# And some helpers
def loglikelihood(*args, **kwargs):
    return log(likelihood(*args, **kwargs))
def nll(*args, **kwargs):
    return -1*loglikelihood(*args, **kwargs)

def normalize(func, **kwargs):
    def f(d):
        return likelihood(d, **kwargs)
    func.norm = 1.0
    func_int, err = integrate.quad(f, -np.inf, np.inf) 
    likelihood.norm = 1.0 / func_int


def get_interval(func, percentage, start, stop, num_points=2000):
    points = np.linspace(start, stop, num_points)
    pair_list = zip(points, map(func, points))

    # ordering rule is maximum likelihood
    # sort by descending in likelihood
    pair_list = sorted(pair_list, key=lambda pair: pair[1], reverse=True)

    delta = (stop-start) / num_points

    total_likelihood = 0.0
    accepted_point_list = []
    for pair in pair_list:
        accepted_point_list.append(pair[0])
        total_likelihood += pair[1]*delta
        #print "Point: ", pair[0], " Likelihood: ", pair[1], " Probability: ", pair[1]*delta, " Total Likelihood", total_likelihood
        if total_likelihood >= percentage: break
        
    interval = (min(accepted_point_list), max(accepted_point_list))
    return interval


# Do the Neyman Construction of a single variable
def get_neyman(func, percentage, param, data):

    interval_list = []

    for param_point in param.linspace():

        def f(d):
            return likelihood(d, **{param.name:param_point})
        likelihood.norm = 1.0
        func_int, err = integrate.quad(f, -np.inf, np.inf) 
        likelihood.norm = 1.0 / func_int
        #f = lambda d: 
        likelihood.normalization = 1.0 / integrate.quad(f, -np.inf, np.inf )[0]
        interval = get_interval(f, percentage, data.min, data.max, data.num_points)
        interval_list.append((param_point, interval))

    return interval_list


def invert_neyman(interval_list, data_point):
    """ Invert neyman to get confidence interval

    The Neyman list looks like: [ (mu, (d0, d1)), ...
    """
    
    mu_list = []

    for item in interval_list:
        (mu, (d0, d1)) = item
        if d0 <= data_point and data_point <= d1:
            mu_list.append(mu)
        pass

    return (min(mu_list), max(mu_list))



def make_plot(mu):
    x = np.linspace(0, 20, 1000) # 100 linearly spaced numbers
    y = [likelihood(d, mu) for d in x]
    plt.plot(x,y)    
    x1,x2,y1,y2 = plt.axis()
    def f(d):
        return likelihood(d, **{"mu":mu})
    likelihood.norm = 1.0
    func_int, err = integrate.quad(f, -np.inf, np.inf) 
    likelihood.norm = 1.0 / func_int
    #likelihood.norm = 1.0
    #likelihood_int, err = integrate.quad(f, -np.inf, np.inf) 
    #likelihood.norm = 1.0 / likelihood_int
    interval = get_interval(f, .68, start=0.0, stop=20.0)
    #print "Interval: ", interval
    plt.vlines(interval[0], y1, y2)
    plt.vlines(interval[1], y1, y2)
    plt.xlabel('x')
    plt.ylabel('likelihood(x)')
    plt.savefig("plot.pdf")



def main():

    # Plot the likelihood as a function of data:
    
    data_var = variable("data", 0, 10, 100)
    mu_var = variable("mu", 2.0, 8.0, 100)

    data_meas = 5

    make_plot(mu=7.39698)

    # clear the current figure
    plt.clf()

    neyman = get_neyman(likelihood, 0.68, mu_var, data_var)
    for pair in neyman:
        #print pair
        (mu, x0, x1) = pair[0], pair[1][0], pair[1][1]
        plt.hlines(mu, x0, x1)

    #print "inverted Neyman 3.0: ", invert_neyman(neyman, 3.0)
    print "inverted Neyman 4.5: ", invert_neyman(neyman, 4.5)
    print "inverted Neyman 5.0: ", invert_neyman(neyman, 5.0)
    print "inverted Neyman 5.5: ", invert_neyman(neyman, 5.5)
    print "inverted Neyman 6.0: ", invert_neyman(neyman, 6.0)
    print "inverted Neyman 8.0: ", invert_neyman(neyman, 8.0)

    plt.xlabel('x')
    plt.ylabel('mu')
    plt.savefig("neyman.pdf")


    # Build the confidence range for the data

    # Build the Neyman Construction over the variable mu0




if __name__=="__main__":
    main()
