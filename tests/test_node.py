
from math import sqrt
from scipy.stats import norm

from variable import *
from node import *
from pdf import pdf

import unittest

def gauss(x, mu, sigma):
    return norm.pdf(x, mu, sigma)

def make_node():
    
    (x0, mu0, sigma0) = make_variables("x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3")
    mass0 = node("mass0", gauss, {'x':x0, 'mu':mu0, 'sigma':sigma0})
    
    (x1, mu1, sigma1) = make_variables("x1[1.2,-5,5], mu1[1,-5,5], sigma1[2,0,3")
    mass1 = node("mass1", gauss, {'x':x1, 'mu':mu1, 'sigma':sigma1})
    
    print "x0 val: ", x0.val, " x1 val: ", x1.val
    
    print "mass0: ", mass0.getVal()
    print "mass1: ", mass1.getVal()
    
    inv_mass = node("inv_mass", invariant_mass, {"a": mass0, "b":mass1})
    
    return inv_mass


def invariant_mass(a, b):
    return sqrt(a*a + b*b)


class TestNode(unittest.TestCase):


    def test_create_vars(self):
        for var in make_variables("x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3], fish"):
            print var, var.name, var.val, var.min, var.max

    def test_node(self):
        inv_mass = make_node()
        print inv_mass.getVal()
        print inv_mass.var('x0').val
    

if __name__ == "__main__":
    unittest.main()
