
from node import node, CannotAnalyticallyIntegrate

from scipy.stats import poisson
from scipy.stats import norm

import sympy

"""
Create common pdf functions that work well with
the PyLikelihood node structure

"""

def _gaussian_distribution(x, mu, sigma):
    return norm.pdf(x, mu, sigma)

def _poisson_distribution(x, mu):
    return poisson.pmf(math.floor(x), mu)


def _gaussian_integral_x(mu, sigma):
    return sympy.N(sympy.pi)*abs(sigma.val)

class gaussian(node):
    """ A Gaussian Probability Density Function

    """

    def __init__(self, name, children):
        print "Initializing Gaussian"
        node.__init__(self, name, _gaussian_distribution, children)
    
    def _analytical_integral(self, *vars_to_integrate):
        """ Implement analytical integral for the gaussian, if possible.

        We can analytically integrate over x and mu, so we
        return this value if we're asked for it.
        Otherwise, we return a numeric integral.
        """
        
        # If we're integrating over x
        if len(vars_to_integrate)==1 and vars_to_integrate[0] == self._children['x']:
            print "Integrating over x"
            if self._children['x'].__class__.__name__ == 'variable':
                print "'x' is a variable"
                return _gaussian_integral_x(self._children['mu'], self._children['sigma'])
        else:
            print "Not integrating over x:"
            print vars_to_integrate
            print self._children['x']

        # If we're integrating over mu
        if vars_to_integrate == [self._children['mu']]:
            print "Integrating over 'mu'"
            if self._children['x'].__class__.__name__ == 'variable':
                return _gaussian_integral_mu(self._children['x'], self._children['sigma'])

        raise CannotAnalyticallyIntegrate

        print "Failed to perform analytic integral"
        raise Exception("AnalyticINtegral")
