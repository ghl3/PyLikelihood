
from node import make_variables
from functions import gaussian

import unittest


class TestGauss(unittest.TestCase):

    def test_gauss_node(self):

        (mass, mu, sigma) = make_variables("mass[10,0,20], mu[0], sigma[1,0,3]")        
        model = gaussian("model_of_mass", [mass, mu, sigma])

    def test_gauss_integral(self):
        
        (mass, mu, sigma) = make_variables("mass[10,0,20], mu[10], sigma[1,0,3]")        
        model = gaussian("model_of_mass", [mass, mu, sigma])

        print model._children

        vars_to_integrate = [mass]
        
        if not vars_to_integrate == [model._children['x']]:
            print model._children['x']
            raise Exception("Integrating variable not in children")

        if not model._children['x'].__class__.__name__ == 'variable':
            print model._children['x'].__class__.__name__
            raise Exception("Integrating variable that is not a 'variable'")

        mass_int = model.integral(mass)
        print "Integral: ", mass_int
        self.assertAlmostEqual(mass_int, 1.0, 6)

