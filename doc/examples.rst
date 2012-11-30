
.. _contents:

.. _library-intro:

Examples
========


Create variables

Variables can be created in groups, and they can be initialized with inital values, ranges, or both.

.. code-block:: python

   (x0, mu0, sigma0) = make_variables("x0, mu0[0], sigma0[0,3], tau[5, 0, 5]")

Create a 'node' out of these variables.  A node is a class that takes a function as well as a set of variables and becomes interpreted as a function that maps the input function over the values of the supplied variables.  It will query their current values and evaluate using those values.


.. code-block:: python

   def gauss(x, mu, sigma):
      return norm.pdf(x, mu, sigma)

   # Variables can be supplied using a dict mapping
   # the function args to the variables  
   (x0, mu0, sigma0) = make_variables("x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3]")
   mass0 = node("mass0", gauss, {'x':x0, 'mu':mu0, 'sigma':sigma0})
    
   # Or they can be supplied as a sequential list,
   # where they are mapped based on order.
   (x1, mu1, sigma1) = make_variables("x1[1.2,-5,5], mu1[1,-5,5], sigma1[2,0,3]")
   mass1 = node("mass1", gauss, [x1, mu1, sigma1])


Nodes can be combined using standard arithmetic expressions.

.. code-block:: python

   mass_product = mass0*mass1


And nodes can be further combined using other functions.


.. code-block:: python

   def invariant_mass(a, b):
      return sqrt(a*a + b*b)

   inv_mass = node("inv_mass", invariant_mass, [mass0, mass1])

