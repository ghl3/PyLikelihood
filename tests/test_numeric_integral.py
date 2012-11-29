
from math import sqrt
from scipy.stats import norm

from node import *
from pdf import pdf

import unittest


def gauss(x, mu, sigma):
    return norm.pdf(x, mu, sigma)

def invariant_mass(a, b):
    return sqrt(a*a + b*b)


class TestIntegral(unittest.TestCase):
    
    
    def setUp(self):

        (x0, mu0, sigma0) = make_variables("x0[.2,-6,6], mu0[0,-5,5], sigma0[1,0,3")
        mass0 = node("mass0", gauss, {'x':x0, 'mu':mu0, 'sigma':sigma0})
    
        (x1, mu1, sigma1) = make_variables("x1[1.2,-6,6], mu1[1,-5,5], sigma1[2,0,3")
        mass1 = node("mass1", gauss, {'x':x1, 'mu':mu1, 'sigma':sigma1})
    
        inv_mass = node("inv_mass", invariant_mass, {"a": mass0, "b":mass1})

        self.vars = (x0, mu0, sigma0, x1, mu1, sigma1)
        self.nodes = (inv_mass, mass0, mass1)


    def test_1d_int(self):
        
        (x0, mu0, sigma0, x1, mu1, sigma1) = self.vars
        (inv_mass, mass0, mass1) = self.nodes

        mu0.val = 0.0
        mu0.sigma = 1.0
        integral = mass0.integral(x0)
        self.assertAlmostEqual(integral, 1.0, 4)
        print "Integral over x0: ", integral

        
    def test_2d_int(self):
        
        (x0, mu0, sigma0, x1, mu1, sigma1) = self.vars
        (inv_mass, mass0, mass1) = self.nodes

        mu0.val = 0.0
        mu0.sigma = 1.0
        mu1.val = 0.0
        mu1.sigma = 1.0

        prod = mass0*mass1

        #integral = inv_mass.integral(x0, x1)
        integral = prod.integral(x0, x1)
        self.assertAlmostEqual(integral, 1.0, 1)
        print "Integral over x0 and x1: ", integral


    def test_3d_int(self):
        
        (x0, mu0, sigma0, x1, mu1, sigma1) = self.vars
        (inv_mass, mass0, mass1) = self.nodes

        (x2, mu2, sigma2) = make_variables("x2[2.2,-6,6], mu2[2,-8,8], sigma2[2,0,3]")
        mass2 = node("mass2", gauss, {'x':x2, 'mu':mu2, 'sigma':sigma2})

        prod = mass0*mass1*mass2

        #integral = inv_mass.integral(x0, x1)
        integral = prod.integral(x0, x1, x2)
        self.assertAlmostEqual(integral, 1.0, 1)
        print "Integral over x0, x1, and x2: ", integral
