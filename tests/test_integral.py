
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

        (x0, mu0, sigma0) = make_variables("x0[.2,-6,6], mu0[0,-5,5], sigma0[1,0,3]")
        mass0 = node("mass0", gauss, {'x':x0, 'mu':mu0, 'sigma':sigma0})
    
        (x1, mu1, sigma1) = make_variables("x1[1.2,-6,6], mu1[1,-5,5], sigma1[2,0,3]")
        mass1 = node("mass1", gauss, {'x':x1, 'mu':mu1, 'sigma':sigma1})
    
        inv_mass = node("inv_mass", invariant_mass, {"a": mass0, "b":mass1})

        self.vars = (x0, mu0, sigma0, x1, mu1, sigma1)
        self.nodes = (inv_mass, mass0, mass1)


    def test_non_dependent_int(self):
        
        (x0, mu0, sigma0, x1, mu1, sigma1) = self.vars
        (inv_mass, mass0, mass1) = self.nodes

        desired_val = mass0()*(x1.max - x1.min)
        integral = mass0.integral(x1)
        self.assertAlmostEqual(desired_val, integral, 4)


    def test_semi_dependent_int(self):
        
        (x0, mu0, sigma0, x1, mu1, sigma1) = self.vars
        (inv_mass, mass0, mass1) = self.nodes

        desired_val = (x1.max - x1.min)
        integral = mass0.integral(x0, x1)
        self.assertAlmostEqual(desired_val, integral, 4)


    def test_product_int(self):
        
        (x0, mu0, sigma0, x1, mu1, sigma1) = self.vars
        (inv_mass, mass0, mass1) = self.nodes

        prod = mass0*mass1

        desired_val = mass1()
        int0 = prod.integral(x0)
        int_mass0 = mass0.integral(x0)
        val_mass1 = mass1()
        print "Desired Val: %s int: %s mass0 int: %s mass1 val: %s" % (desired_val, int0, int_mass0, val_mass1) 
        self.assertAlmostEqual(int0, desired_val, 3)

        desired_val = mass0()
        int1 = prod.integral(x1)
        int_mass1 = mass1.integral(x1)
        val_mass0 = mass0()
        print "Desired Val: %s int: %s mass1 int: %s mass0 val: %s" % (desired_val, int1, int_mass1, val_mass0) 
        self.assertAlmostEqual(int1, desired_val, 3)

        int1 = prod.integral(x1)*mass0()
        integral = prod.integral(x0, x1)
        print "int0: %s, int1: %s, total: %s" % (int0, int1, integral)
        self.assertAlmostEqual(integral, 1.0, 1)


