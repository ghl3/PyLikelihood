#!/usr/bin/env python

from __future__ import division

import inspect

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import poisson
from scipy.stats import norm
from scipy import integrate

from math import log

class variable(object):
    """ A class to store a name, value, and range

    """

    def __init__(self, name, var_min=-np.inf, var_max=np.inf, num_points=1000):
        self.name = name
        self.min = var_min
        self.max = var_max
        self.num_points = num_points
        self.val = 0.0

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

    def __init__(self, pdf, data=None, params=None):
        self.pdf = pdf

        # Use inspection to find all arguments of the pdf
        func_spec = inspect.getargspec(self.pdf)
        (all_arguments, all_defaults) = (func_spec.args, func_spec.defaults)

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
        setattr(self, data.name, data)

        # Get the list of parameters,
        # create variables based on them,
        # and store that list of variables
        self.param_list = []
        param_list = [arg for arg in func_spec.args if arg != self.data.name]
        print "Param List: ", param_list 
        for param in param_list:

            # Check if the parameter matches one in the
            # supplied list of parameter 'variables'
            matching_var_list = [var for var in params if var.name==param]
            
            param_var = None

            if len(matching_var_list)==0:
                print "Creating new param with name: ", param
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
            self.param_list.append(param_var)
            # And create an attribute for easy access
            print param_var, param_var.name
            setattr(self, param_var.name, param_var)
        
        self.norm = 1.0
        self.normalization_cache = {}


    def state(self):
        """ Return a dict with the current state of data and parameters

        """
        current_state = {}
        #current_state[self.data.name] = self.data.val
        for param in self.param_list:
            current_state[param.name] = param.val
        return current_state


    def _eval_raw(self, x):
        """ Get the current state of the likelihood
        without any normalization

        """
        current_state = self.state()
        return self.pdf(x, **current_state)
        

    def eval(self, x):
        """ Val of pdf based on the current state,
        Evaluated on the given data point
        This includes normalization, which is cached

        """
        self.normalize()
        return self._eval_raw(x)*self.norm 


    def normalize(self):
        """ Integrate over the data
        at the current parameter point
        """

        # Check if the normalization has been cached
        # If so, return that cache
        param_state = hash(frozenset(self.state().items()))
        try:
            #print "Found Norm in Cache"
            return self.normalization_cache[param_state]
        except KeyError:
            pass

        # If not, integrate over the data, invert it, 
        # and store the cache
        data_min, data_max = (self.data.min, self.data.max)
        print "Integrating: ", self.state()
        integral, err = integrate.quad(self._eval_raw, data_min, data_max) 
        self.norm = 1.0 / integral
        self.normalization_cache[param_state] = self.norm
        return self.norm
    

    # Make the class callable:
    def __call__(self, *args, **kwargs):
        return self.eval(*args, **kwargs)


    # Log Likelihood
    def loglikelihood(self, *args, **kwargs):
        return log(self.eval(*args, **kwargs))


    # Negative Log Likelihood
    def nll(self, *args, **kwargs):
        return -1*self.loglikelihood(*args, **kwargs)

    
    def get_interval(self, percentage):
        """ Get an inverval over the data parameter 
        which contains 'percentage' of the probability
        
        """
        points = self.data.linspace()
        pair_list = zip(points, map(self.eval, points))

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

        param_var = getattr(self, param)
        
        interval_list = []
        for param_point in param_var.linspace():
            param_var.val = param_point
            interval = self.get_interval(percentage)
            interval_list.append( (param_point, interval) )

        return interval_list


    def make_plot(self, interval=None):
        """ Plot the likelihood over data
        """
        x = self.data.linspace()
        y = map(self.eval, x)
        plt.plot(x,y)
        plt.xlabel('x')
        plt.ylabel('likelihood(x)')
        x1,x2,y1,y2 = plt.axis()

        if interval != None:
            interval = self.get_interval(interval)
            plt.vlines(interval[0], y1, y2)
            plt.vlines(interval[1], y1, y2)

        return


    def invert_neyman(self, data_point, neyman=None, percentage=None, param=None):
        """ Invert neyman to get confidence interval
        
        The Neyman list looks like: [ (mu, (d0, d1)), ...
        """
        
        mu_list = []
        
        if neyman==None:
            neyman = self.get_neyman(percentage, param)

        for item in neyman:
            (mu, (d0, d1)) = item
            if d0 <= data_point and data_point <= d1:
                mu_list.append(mu)
            pass

        return (min(mu_list), max(mu_list))


def main():

    # Annoying boilerplate
    def pois(x, mu):
        return poisson.pmf(x, mu)
    def gauss(x, mu, sigma):
        return norm.pdf(x, mu, sigma)

    d = variable("d", 0, 10, 100)
    mu = variable("mu", 0, 10, 100)
    mu0 = variable("mu0", 0, 10, 100)

    def pdf(d, mu=5.0, mu0=5.0, sigma=2.0):
        return gauss(d, mu, 2.0)*gauss(mu0, mu, sigma)

    # Create the likelihood
    model = likelihood(pdf, data=d, params=[mu, mu0])

    print model

    from pprint import pprint
    pprint (vars(model))
    print model.mu
    print model.mu0

    model.mu.val = 5.0
    model.mu0.val = 5.0
    model.sigma.val = 1.0

    print model.norm
    model.eval(5)
    print model.norm

    for point in range(0, 10):
        print point, ": ", model(point), ": ", model.eval(point)

    print model.get_interval(0.68)

    model.make_plot(interval=0.68)
    plt.savefig("plot.pdf")

    # clear the current figure and
    # plot the neyman interval
    plt.clf()
    neyman = model.get_neyman(0.68, "mu")
    for pair in neyman:
        #print pair
        (mu, x0, x1) = pair[0], pair[1][0], pair[1][1]
        plt.hlines(mu, x0, x1)
    plt.xlabel('x')
    plt.ylabel('mu')
    plt.savefig("neyman.pdf")

    print "inverted Neyman 4.5: ", model.invert_neyman(4.5, neyman)
    print "inverted Neyman 5.0: ", model.invert_neyman(5.0, neyman)
    print "inverted Neyman 5.5: ", model.invert_neyman(5.5, neyman)
    print "inverted Neyman 6.0: ", model.invert_neyman(6.0, neyman)
    print "inverted Neyman 8.0: ", model.invert_neyman(8.0, neyman)


if __name__=="__main__":
    main()
