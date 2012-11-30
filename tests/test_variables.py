

from node import *

import numpy as np
import unittest


class TestVariables(unittest.TestCase):


    def test_create_vars(self):

        (x0, mu0, sigma0, fish, frank) = make_variables("x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3] fish frank")
        
        self.assertEqual(x0.name, 'x0')
        self.assertEqual(x0.min, -5)
        self.assertEqual(x0.max, 5)
        self.assertEqual(x0.val, .2)

        self.assertEqual(fish.name, 'fish')
        self.assertEqual(fish.min, -np.inf)
        self.assertEqual(fish.max, np.inf)

        self.assertEqual(frank.name, 'frank')
        self.assertEqual(frank.min, -np.inf)
        self.assertEqual(frank.max, np.inf)



    def test_bad_vars(self):

        with self.assertRaises(Exception):
            (x0, mu0, sigma0) = make_variables("x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3")

