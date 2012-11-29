

from math import sqrt, exp, factorial
from scipy.stats import norm

from node import *
from pdf import pdf

import unittest


def pois(k, lam):
    return lam**k*exp(-1*lam)/factorial(k)


class TestExample1(unittest.TestCase):
    
    
    def test_example(self):
        (d, s, b) = make_variables("d[1, 0, 15], s[1,0,15], b[1, 0, 15]")
        chan = node("chan", pois, [d, N_exp])
        
        (x0, mu0, sigma0) = make_variables("x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3]")
        mass0 = node("mass0", gauss, {'x':x0, 'mu':mu0, 'sigma':sigma0})
        
        (x1, mu1, sigma1) = make_variables("x1[1.2,-5,5], mu1[1,-5,5], sigma1[2,0,3]")
        mass1 = node("mass1", gauss, {'x':x1, 'mu':mu1, 'sigma':sigma1})
        
        print "x0 val: ", x0.val, " x1 val: ", x1.val
        
        print "mass0: ", mass0.getVal()
        print "mass1: ", mass1.getVal()
        
        inv_mass = node("inv_mass", invariant_mass, {"a": mass0, "b":mass1})
        
        return inv_mass
