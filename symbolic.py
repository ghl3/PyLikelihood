

import sympy


class poisson(sympy.Function):

     nargs = 1

      @classmethod
      def eval(cls, k, lam):
          if k >= 0:
              return lam**k * sympy.exp(-1*lam) / sympy.factorial(k)
          else:
              return 0

      def fdiff(self, argindex=1):
          """
          Returns the first derivative of this function.
          """
          if argindex == 1:
              return 
          else:
              raise ArgumentIndexError(self, argindex)

    return 


def main():
    """ Testing SymPy functionality

    Hoping to add symbolic functions and integrals

    """


    
