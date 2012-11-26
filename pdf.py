
import inspect
import itertools

from node import node
from variable import variable

from scipy import integrate
from scipy import optimize

from math import log
from math import ceil

import logging
logging.basicConfig()


class pdf(object):
    """ A function which evaluates a pdf on data

    The pdf class is instantiated by a function, and it
    takes care of ensuring that the function is properly
    normalized by integrating over that function.

    It is specialized for the case where that function
    is a node tree and uses that node's specific integration
    capabilities.
    
    The first argument of the function is considered to be the "data"
    (it can be a list, of course) and the rest of the function
    is assumed to be parameters.

    The pdf object contains an internal cache to speed up
    its process of normalizing itself, and also tries to use
    either symbolic or node structure of its dependencies
    for normalization (if it can)
    """

    # Create an internal logger
    logging = logging.getLogger("pdf")

    def __init__(self, func, data):
        """ Create a pdf from a function

        The function must be an instance of a 'node' class

        In the construction, one must specify which variables
        are to be interpreted as data.  This can be done by
        supplying a list of variables or a string list of their
        names.

        """

        self.name = "pdf " + func.name
        self._func = func
        self._data = []

        # Normalization
        self.norm = 1.0
        self.normalization_cache = {}
        self.minimization_cache = {}

        # Always insist that we have data_vars supplied
        if data==None or len(data)==0:
            print "Error: No Data Vars Supplied"
            raise Exception()

        # Handle the special 'node' case
        if func.__class__.__name__ != "node":
            print "(For now) All pdf's must be made from nodes"
            raise Exception()
            
        # Get all dependencies
        all_arguments = func.dependent_vars()
        for var in data:
            if var.__class__.__name__ == "variable":
                self._data.append(var)
            else:
                var = func.var(var)
                self._data.append(var)
        self._params = [var for var in all_arguments
                        if var not in self._data]
        pass


    def _eval_raw(self):
        """ Get the current value of the pdf w/o normalization

        This will use all current values of the dependent
        variables and the data.
        """

        pdf_val = self._func()
        if pdf_val < 0:
            print "Error: Pdf is 0 at state: ", current_state
            raise Exception("PdfVal")
        return pdf_val


    def eval(self, **kwargs):
        """ Get the current value of the pdf

        This value is determined by evaluating the
        pdf function on the current state of parameters
        and data, including normalization over all data.
        
        The values of parameters or data can be set
        by specifying their values as keyword args.
        """

        self.set_state(**kwargs)

        self.normalize()

        likelihood_val = self._eval_raw()*self.norm 
        if likelihood_val < 0:
            print "Error: Pdf evaluates to <0 at state: ", current_state
            raise Exception("Pdf val is negative")

        return likelihood_val


    # Make the class callable:
    def __call__(self, *args, **kwargs):
        return self.eval(*args, **kwargs)

    
    def var(self, var_name):
        """ Return the dependent var by name

        """

        for var in itertools.chain(self._params, self._data):
            if var.name == var_name:
                return var

        # If we get here, then we didn't find the variable
        print "Error: Didn't find variable: %s in pdf: %s" % (var_name, self.name)
        raise KeyError("Variable %s not found" % var_name)


    def total_state(self):
        """ Return a dict with the current state of all variables

        This includes both parameters and data variables
        """

        current_state = {}
        for param in itertools.chain(self._params, self._data):
            current_state[param.name] = self.var(param.name).val
        #for param in self._data:
        #    current_state[param.name] = self.var(param.name).val
        return current_state


    def param_state(self):
        """ Return a dict with the current state only the parameters

        """
        current_state = {}
        for param in self._params:
            current_state[param.name] = self.var(param.name).val
        return current_state


    def set_state(self, **kwargs):
        """ Set the state based on the supplied keyword args

        It will throw an exception if any unknown variables
        are attempted to be set
        """

        for (arg, val) in kwargs.iteritems():
            self.var(arg).val = val


    def set_data(self, data):
        """ Set the value of the data parameter(s)

        This takes either a single val, a list, 
        or a dictionary
        """

        # Dictioanry
        try:
            for var, val in data.iteritems():
                var.val = val
            return
        except AttributeError:
            pass

        # List
        try:
            for var, val in zip(self._data, data):
                var.val = val
            return
        except TypeError:
            pass
        
        # Single entry
        if len(self._data)==1:
            self._data[0].val = data
            return

        # If we get here, we failed pretty hard
        print "Error: Cannot set data based on input: ", data
        raise Exception()


    def normalize(self):
        """ Determine the integral of the pdf over data

        This integral is evaluated over all data variables
        and is evaluated at the current state of all parameters.
        Normalization integrals are cached.
        """

        # Check if the normalization has been cached
        # If so, return that cache
        norm_key = str(self.param_state().items())
        try:
            norm = self.normalization_cache[norm_key]
            self.norm = norm
            print "Using Cached Normalization for %s: %s" % (self.name, norm)
            return norm
        except KeyError:
            pass

        # We have to normalize, but we should be sure
        # to restore the state after, so we aren't effected
        # by the random data points used to evaluate the integral
        data_before = {var.name : var.val for var in self._data}
        #for var in self._data:
        #    data_before[var.name] = var.val

        # To normalize, we integrate the 'raw' function
        # over data

        if len(self._data)==1:
            def func_for_norm(data_val):
                self._data[0].val = data_val
                return self._eval_raw()
            data_min, data_max = (self._data[0].min, self._data[0].max)
            integral, err = integrate.quad(func_for_norm, data_min, data_max) 
            self.norm = 1.0 / integral

        elif len(self._data)==2:
            def func_for_norm(data1_val, data0_val):
                self._data[0].val = data0_val
                self._data[1].val = data1_val
                return self._eval_raw()
            data0_min, data0_max = (self._data[0].min, self._data[0].max)
            data1_min, data1_max = (self._data[1].min, self._data[1].max)
            integral, err = integrate.dblquad(func_for_norm, data0_min, data0_max,
                                              lambda x: data1_min, lambda x: data1_max)
            self.norm = 1.0 / integral
        
        else:
            raise Exception("Data of dim > 2 not currently handled")

        if self.norm <= 0:
            print "Error: Normalization is <= 0"
            raise Exception("BadNormalization")

        #print "Found integral: ", integral, " Normalization: ", self.norm
        self.normalization_cache[norm_key] = self.norm

        self.logging.debug("Got Norm From Integral: %s from state: %s" % \
                               (self.norm, str(self.param_state())) )

        # Restore the data values
        for var in self._data:
            var.val = data_before[var.name]

        return self.norm


    def fitTo(self, data, params_to_fit, **kwargs):
        """ Minmize the supplied parameters based on the nll

        Set the values of the minimized parameters in the
        likelihood's 'state'.  Use any keyword arguments as
        initial values to parameters, or for any other
        (non-minimized) parameters in the model

        TO DO: Split the kwargs into args for the nll
        and args for optimize.minimize, and throw
        exceptions for all others
        (Should also warn about argument clashes with
        optimize when the likelihood function is initialized...)

        """

        # Log
        self.logging.debug( "Minimizing: " + str(params_to_fit) 
                            + " on state: " + str(self.total_state()))

        # Minimize the supplied params
        if len(params_to_fit)==0:
            print "Error: No Paramaterize to Minimize"
            raise Exception("FitTo")

        # Set the value of the data to the supplied data
        self.set_data(data)

        # NOT YET IMPLEMENTED
        # Create a key based on the values of the params to not minimize
        # and the list of params to minimize (possibly overkill, but whatevs)
        # as well as the data being minimized
        constant_params = [item for item in self.param_state().items() 
                           if item[0] not in params_to_fit]

        cache_key = (tuple(self._data), frozenset(constant_params), frozenset(params_to_fit))
        if cache_key in self.minimization_cache:
            state = self.minimization_cache[cache_key] 
            print "Using fitTo Cache: ", state
            self.set_state(**state)
            return

        #current_state = self.total_state()
        self.normalize()

        # Create the function for minimization
        def nnl_for_min(param_value_list):
            """ Create the wrapper function for scipy.optimize

            """

            for (param, val) in zip(params_to_fit, param_value_list):
                self.var(param).val = val

            # nll without normalization
            return -1*log(self._eval_raw())


        # Get the initial guess
        guess = [self.var(param).val for param in params_to_fit]

        # Run the minimization
        bounds = []
        for param in params_to_fit:
            param_min = self.var(param).min
            param_max = self.var(param).max
            #param_min = getattr(self, param+"_var").min 
            #param_max = getattr(self, param+"_var").max
            bounds.append( (param_min, param_max) )
        #res = optimize.minimize(nnl_for_min, guess, method="TNC", 
        #                        bounds=bounds, tol=0.00000001)
        res = optimize.minimize(nnl_for_min, guess, method="BFGS")


        self.logging.debug("Successfully Minimized the function: " + str(res))
        
        # Take the result and set all parameters
        def set_to_minimum(res):
            min_values = res.x
            for (param, val) in zip(params_to_fit, min_values):
                param_min = self.var(param).min
                param_max = self.var(param).max
                #param_min = getattr(self, param+"_var").min 
                #param_max = getattr(self, param+"_var").max
                if val < param_min : val = param_min
                if val > param_max : val = param_max
                self.logging.debug("Minimized value of %s : %s" % (param, val))
                #setattr(self, param, val)
                self.var(param).val = val

        set_to_minimum(res)

        # Cache the result
        self.minimization_cache[cache_key] = self.total_state()

        return

