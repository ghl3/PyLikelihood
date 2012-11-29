
from math import sqrt
from scipy.stats import norm

#from variable import *
from node import *
from pdf import pdf

import unittest


class TestArithmetic(unittest.TestCase):

    def setUp(self):
        pass

    def test_variable_sum(self):
        (a, b, c, d) = make_variables("a[1], b[2], c[3], d[4]")
        a_plus_b = a+b
        print a_plus_b
        self.assertEqual(a_plus_b(),  a() + b())
        self.assertEqual(a_plus_b(),  3)

    def test_add(self):
        (a, b, c, d) = make_variables("a[1], b[2], c[3], d[4]")
        a_plus_b = sum_node("a_plus_b", a, b)
        print a_plus_b
        self.assertEqual(a_plus_b(),  a() + b())
        self.assertEqual(a_plus_b(),  3)

    def test_subtract(self):
        (a, b, c, d) = make_variables("a[1], b[2], c[3], d[4]")
        a_minus_b = diff_node("a_minus_b", a, b)
        print a_minus_b
        self.assertEqual(a_minus_b(), a() - b())
        self.assertEqual(a_minus_b(), -1)

    #c_plus_d = c*d
    #c_div_d = c/d

if __name__ == "__main__":
    unittest.main()

