

from math import sqrt, exp, factorial
from scipy.stats import norm

from node import *
from pdf import pdf

import unittest


def pois(k, lam):
    return lam**k*exp(-1*lam)/factorial(k)

def gauss(x, mu, sigma):
    return norm.pdf(x, mu, sigma)


class TestExample1(unittest.TestCase):
    
    
    def test_example(self):
        (d, s, b) = make_variables("d[1, 0, 15], s[1,0,15], b[1, 0, 15]")
        chan = node("chan", pois, [d, s+b])
        
        (b0, sigma_b) = make_variables("b0[5,0,10], sigma_b[1,0,3]")
        b_constraint = node("b_constraint", gauss, [b0, b, sigma_b])
        
        model = chan*b_constraint
        
        print "model val: ", model()

