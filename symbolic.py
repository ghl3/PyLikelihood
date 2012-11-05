
from sympy import *

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

    sym_norm = integrate(poisson(d, 5+3)*gauss(5, 5, 1.0), (d, 0, 100))
    print sym_norm

    norm = N(integrate(poisson(d, 5+3)*gauss(5, 5, 1.0), (d, 0, 100)))
    print norm



if __name__ == "__main__":
     main()
