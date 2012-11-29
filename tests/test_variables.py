

from node import *

import unittest

class TestVariables(unittest.TestCase):
    

    def test_create_vars(self):

        with self.assertRaises(Exception):
            (x0, mu0, sigma0) = make_variables("x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3")
