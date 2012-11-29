

from pdf import pdf

from node import node
from variable import make_variables
from test_node import gauss, make_node, invariant_mass

import logging


import unittest


class TestPdf(unittest.TestCase):

    def test_vals(self):
        
        # Test a gaussian pdf
        (x, mu, sigma) = make_variables("x[.2,-5,5], mu[0,-5,5], sigma[1,0,3")
        mass = node("mass", gauss, {'x':x, 'mu':mu, 'sigma':sigma})
        my_pdf = pdf(mass, data=['x'])
        
        # Evaluate the pdf    
        mu.val = 0    
        val = my_pdf()
        print "x=%s, mu=%s, sigma=%s, val=%s" % (x.val, mu.val, sigma.val, val)
        self.assertTrue( abs(val - 0.39104269) < .01 )
        
        mu.val = 2    
        val = my_pdf()
        print "x=%s, mu=%s, sigma=%s, val=%s" % (x.val, mu.val, sigma.val, val)
        self.assertTrue( abs(val - 0.07895016) < .01 )

        x.val = 1    
        val = my_pdf()
        print "x=%s, mu=%s, sigma=%s, val=%s" % (x.val, mu.val, sigma.val, val)
        self.assertTrue( abs(val - 0.24197072) < .01 )


    def test_pdf(self):

        # Make a simple pdf for integration and normalization testing
        (x0, mu0, sigma0) = make_variables("x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3")
        mass0 = node("mass0", gauss, {'x':x0, 'mu':mu0, 'sigma':sigma0})
        my_pdf = pdf(mass0, data=['x0'])

        # Evaluate the pdf
        val = my_pdf()
        print "x=%s, mu=%s, sigma=%s, val=%s" % (x0.val, mu0.val, sigma0.val, val)
        
        val = my_pdf()
        print "x=%s, mu=%s, sigma=%s, val=%s" % (x0.val, mu0.val, sigma0.val, val)
        
        # 1-d Fitting
        my_pdf.logging.setLevel(logging.DEBUG)
        my_pdf.fitTo( data=1, params_to_fit=["mu0"])
        print "Fitted pdf to x=1: mu=%s x=%s" % (mu0.val, x0.val)
        
        # Get a 2-d function and make a pdf
        inv_mass = make_node()
        my_pdf = pdf(inv_mass, data=['x0', 'x1'])

        # By Variable
        x0 = inv_mass.var('x0')
        x1 = inv_mass.var('x1')
        my_pdf = pdf(inv_mass, data=[x0, x1])

        print "Params: ", [param.name for param in my_pdf._params]
        print "Data: ", [param.name for param in my_pdf._data]
        print "x0 in pdf: ", my_pdf.var('x0'), my_pdf.var('x0').val

        # Evaluate the pdf
        val = my_pdf()
        print "x=%s, mu=%s, sigma=%s, val=%s" % (x0.val, mu0.val, sigma0.val, val)

        my_pdf.fitTo( data={x0:1, x1:1.2}, params_to_fit=["mu0"])

if __name__ == "__main__":
    unittest.main()
