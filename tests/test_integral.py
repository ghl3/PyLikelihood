
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

        (x0, mu0, sigma0) = make_variables("x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3")
        mass0 = node("mass0", gauss, {'x':x0, 'mu':mu0, 'sigma':sigma0})
    
        (x1, mu1, sigma1) = make_variables("x1[1.2,-5,5], mu1[1,-5,5], sigma1[2,0,3")
        mass1 = node("mass1", gauss, {'x':x1, 'mu':mu1, 'sigma':sigma1})
    
        inv_mass = node("inv_mass", invariant_mass, {"a": mass0, "b":mass1})

        self.vars = (x0, mu0, sigma0, x1, mu1, sigma1)
        self.node = inv_mass


    def test_1d_int(self):
        
        (x0, mu0, sigma0, x1, mu1, sigma1) = self.vars
        self.node.integral(x0)
        
