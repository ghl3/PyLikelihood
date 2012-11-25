
import inspect
import itertools

from node import node
from variable import variable

'''

THIS IS JUST A TEST, FOR NOW


'''


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
    #logging = logging.getLogger("likelihood")

    def __init__(self, func, data, params=None):

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
        if func.__class__.__name__ == "node":
            
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
            
        # If it is any other type of function
        else:
            print "(For now) All pdf's must be made from nodes"
            raise Exception()
        '''
            func_spec = inspect.getargspec(func)
            (all_arguments, all_defaults) = (func_spec.args, func_spec.defaults)

            self._data = []
            for var in data:
                if var.__class__.__name__ == "variable":
                    if var.name not in all_arguments:
                        print "Error: Supplied data var: ", var.name
                        print " does not match a function arugment for: ", func
                        raise Exception()
                    self._data.append(var)

        if params==None:
            self._params = [variable(arg) for arg in all_arguments[1:]]
        else:
            for var in params:
                if var.__class__.__name__ != "variable":
                    print "Error: Suppied data must be a variable"
                    raise Exception()
            self._params = params

        # Check that all arguments are used and that each
        # argument is used only once
        for arg in all_arguments:
            if arg not in [var.name for var in self._params]:
                if arg not in [var.name for var in self._data]:
                    print "Error: Unhandled argument: ", arg
                    raise Exception()
            if arg in [var.name for var in self._params]:
                if arg in [var.name for var in self._data]:
                    print "Error: Arg is set to be both data and a param ", arg
                    raise Exception()


        '''
        pass


    def _eval_raw(self, data, **kwargs):
        """ Get the current state of the likelihood
        without any normalization

        Any values can be set via kwargs, and the evaluation
        will take place after those kwargs are set
        """

        self.set_state(**kwargs)
        current_state = self.total_state()
        #return self.pdf(self.get_data(), **current_state)
        pdf_val = self.func(data, **current_state)
        if pdf_val < 0:
            print "Error: Pdf is 0 at state: ", current_state
            raise Exception("PdfVal")
        return pdf_val


    def eval(self, data, **kwargs):
        """ Val of pdf based on the current state,
        Evaluated on the given data point
        This includes normalization, which is cached

        The current state can be set via kwargs, and
        the normalization will take place after the
        current state is set.

        """
        self.set_state(**kwargs)

        # We have to normalize, but we should be sure
        # to restore the state after, so we aren't effected
        # by the random data points used to evaluate the integral
        #current_state = self.total_state()
        #current_data = self.get_data()
        self.normalize()

        likelihood_val = self._eval_raw(data)*self.norm 
        if likelihood_val < 0:
            print "Error: Pdf is 0 at state: ", current_state
            raise Exception("PdfVal")
        return likelihood_val


    # Make the class callable:
    def __call__(self, *args, **kwargs):
        return self.eval(*args, **kwargs)

    
    def var(self, var_name):
        for var in itertools.chain(self._params, self._data):
            #print "Checking var: ", var, var.name
            if var.name == var_name:
                return var
        print "Error: Didn't find variable: ", var_name,
        print " in pdf: ", self.name
        raise Exception()


    def normalize(self):
        """ Integrate over the data
        at the current parameter point
        """

        # Check if the normalization has been cached
        # If so, return that cache
        norm_key = str(self.param_state().items())
        try:
            norm = self.normalization_cache[norm_key]
            self.norm = norm
            return norm
        except KeyError:
            pass

        # To normalize, we integrate the 'raw' function
        # over data
        def func_for_norm(data):
            return self._eval_raw(data)

        # If not, integrate over the data, invert it, 
        # and store the cache
        data_min, data_max = (self.data.min, self.data.max)
        #self.logging.debug("Integrating: " + str(self.param_state()))
        #integral, err = integrate.quad(self._eval_raw, data_min, data_max) 
        integral, err = integrate.quad(func_for_norm, data_min, data_max) 
        self.norm = 1.0 / integral

        if self.norm <= 0:
            print "Error: Normalization is <= 0"
            raise Exception("BadNormalization")

        #print "Found integral: ", integral, " Normalization: ", self.norm
        self.normalization_cache[norm_key] = self.norm

        self.logging.debug("Got Norm From Integral: %s from state: %s" % \
                               (self.norm, str(self.param_state())) )

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

        # NOT YET IMPLEMENTED
        # Create a key based on the values of the params to not minimize
        # and the list of params to minimize (possibly overkill, but whatevs)
        # as well as the data being minimized
        constant_params = [item for item in self.param_state().items() 
                           if item[0] not in params_to_fit]

        cache_key = (data, frozenset(constant_params), frozenset(params_to_fit))
        if cache_key in self.minimization_cache:
            state = self.minimization_cache[cache_key] 
            print "Using fitTo Cache: ", state
            self.set_state(**state)
            return

        # Minimize the supplied params
        if len(params_to_fit)==0:
            print "Error: No Paramaterize to Minimize"
            raise Exception("FitTo")
            return

        # Set the value of the data to the supplied data
        self.set_data(data)

        #current_state = self.total_state()
        self.normalize()

        # Create the function for minimization
        def nnl_for_min(param_value_list):
            """ Create the wrapper function for scipy.optimize

            """

            for (param, val) in zip(params_to_fit, param_value_list):
                setattr(self, param, val)

            # nll without normalization
            return -1*log(self._eval_raw())


        def set_to_minimum(res):
            min_values = res.x
            for (param, val) in zip(params_to_fit, min_values):
                param_min = getattr(self, param+"_var").min 
                param_max = getattr(self, param+"_var").max
                if val < param_min : val = param_min
                if val > param_max : val = param_max
                self.logging.debug("Minimized value of %s : %s" % (param, val))
                setattr(self, param, val)

        # Get the initial guess
        guess = [getattr(self, param) for param in params_to_fit]

        # Run the minimization
        bounds = []
        for param in params_to_fit:
            param_min = getattr(self, param+"_var").min 
            param_max = getattr(self, param+"_var").max
            bounds.append( (param_min, param_max) )
        res = optimize.minimize(nnl_for_min, guess, method="TNC", 
                                bounds=bounds, tol=0.000001)

        set_to_minimum(res)

        self.logging.debug("Successfully Minimized the function: " + str(res))

        # Cache the result
        self.minimization_cache[cache_key] = self.total_state()

        return

