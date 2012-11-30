
from math import sqrt
from scipy.stats import norm

from node import *
from pdf import pdf

import unittest


class TestArithmetic(unittest.TestCase):

    def setUp(self):
        pass


    def test_add(self):
        (a, b, c, d) = make_variables("a[1], b[2], c[3], d[4]")
        a_plus_b = sum_node("a_plus_b", a, b)
        print a_plus_b
        self.assertEqual(a_plus_b(),  a() + b())
        self.assertEqual(a_plus_b(),  3)
        a_plus_b = a+b
        print a_plus_b
        self.assertEqual(a_plus_b(),  a() + b())
        self.assertEqual(a_plus_b(),  3)


    def test_subtract(self):
        (a, b, c, d) = make_variables("a[1], b[2], c[3], d[4]")
        a_minus_b = diff_node("a_minus_b", a, b)
        print a_minus_b
        self.assertEqual(a_minus_b(), a() - b())
        self.assertEqual(a_minus_b(), -1)
        a_minus_b = a - b
        print a_minus_b
        self.assertEqual(a_minus_b(), a() - b())
        self.assertEqual(a_minus_b(), -1)


    def test_multiply(self):
        (a, b, c, d) = make_variables("a[1], b[2], c[3], d[4]")
        a_times_b = product_node("a_times_b", a, b)
        print a_times_b
        self.assertEqual(a_times_b(), a() * b())
        self.assertEqual(a_times_b(), 2)
        a_times_b = a*b
        self.assertEqual(a_times_b(), a() * b())
        self.assertEqual(a_times_b(), 2)


    def test_division(self):
        (a, b, c, d) = make_variables("a[1], b[2], c[3], d[4]")
        a_over_b = quotient_node("a_over_b", a, b)
        print a_over_b
        self.assertEqual(a_over_b(), a() / b())
        self.assertEqual(a_over_b(), .5)
        a_over_b = a/b
        print a_over_b
        self.assertEqual(a_over_b(), a() / b())
        self.assertEqual(a_over_b(), .5)


if __name__ == "__main__":
    unittest.main()

