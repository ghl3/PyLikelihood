

from scipy import integrate as sci_integrate
from sympy import *

from Test import pdf, create_model
from time import time

class poisson(Function):

     nargs = 2
     
     @classmethod
     def eval(cls, k, lam):
          if k >= 0:
               return lam**k * exp(-1*lam) / factorial(k)
          else:
               return 0

     '''
     def fdiff(self, argindex=1):
          """
          Returns the first derivative of this function.
          """
          if argindex == 1:
               return 
          else:
               raise ArgumentIndexError(self, argindex)
     '''

class gauss(Function):

     nargs = 3
     
     @classmethod
     def eval(cls, x, mu, sigma):
          return 1/(sigma*sqrt(2*pi)) * exp(-(x-mu)**2 / (2*sigma**2))


          
def main():
    """ Testing SymPy functionality

    Hoping to add symbolic functions and integrals

    """

    n = Symbol('n')
    lam = Symbol('lam')
    
    print poisson(n, lam)
    print N(poisson(1, 1))
    print N(integrate(poisson(n, 1), (n, 0, 10)))

    x, mu, sigma = symbols('x mu sigma')

    print gauss(x, mu, sigma)
    print N(gauss(0, 0, 1))
    print N(integrate(gauss(x, 0, 1), (x, -5, 5)))

    d, s, b, b0, sigma = symbols("d s b b0 sigma")

    on_off = poisson(d, s+b)*gauss(b0, b, sigma)
    print on_off

    print "Eval: "
    print N(on_off, subs={d:1, s:3, b:5, b0:5, sigma:1})
    #print N(on_off.evalf(subs={d:1, s:3, b:5, b0:5, sigma:1})
    return


    sym_norm = integrate(poisson(d, 5+3)*gauss(5, 5, 1.0), (d, 0, 100))
    print sym_norm

    norm = N(integrate(poisson(d, 5+3)*gauss(5, 5, 1.0), (d, 0, 100)))
    print norm


    pdf = create_model()
    data_min, data_max = (pdf.data.min, pdf.data.max)

    # Timing the integral
    int_0_begin = time()
    norm = N(integrate(poisson(d, 5+3)*gauss(5, 5, 1.0), (d, data_min, data_max)))
    int_0_end = time()

    # Integrate the model
    def func_for_norm(data_val):
         pdf.set_data(data_val)
         return pdf._eval_raw()

    int_1_begin = time()
    integral, err = sci_integrate.quad(func_for_norm, data_min, data_max) 
    int_1_end = time()

    print "Symbolic: ", int_0_end-int_0_begin,
    print " Numeric: ", int_1_end-int_1_begin,



if __name__ == "__main__":
     main()
