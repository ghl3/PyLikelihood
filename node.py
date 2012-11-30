
# cause, duh:
from __future__ import division

import inspect

import re
import scipy
import numpy as np


class variable(object):
    """ A class to store a name, value, and range

    """
    # Constructor
    def __init__(self, name, val=None, var_min=-np.inf, var_max=np.inf, num_points=1000):
        self.name = name
        self.min = var_min
        self.max = var_max
        self.num_points = num_points
        self.val = val

    def getVal(self):
        return self.val

    def __call__(self):
        return self.val

    def __add__(self, other):
        """ Return the sum of two variables as a node

        """
        name = self.name + "_plus_" + other.name
        return sum_node(name, self, other)

    def linspace(self):
        return np.linspace(self.min, self.max, self.num_points)


class node(object):
    """ An node in an acyclic graph

    Can be connected to other nodes to form a graph/tree.
    It uses this tree structure to determine its caching
    and its integration

    """

    def __init__(self, name, func, children):

        self.name = name
        self._func = func

        # graph properties
        self._children = {}
        
        # caching properties
        self._var_values = {}
        self._cached_val=None
        self._tolerence = 0.00000001
        self._dirty = True

        func_spec = inspect.getargspec(func)
        (all_arguments, all_defaults) = (func_spec.args, func_spec.defaults)

        # Add the child nodes as children to this node

        # I 'children' is a dict
        if isinstance(children, dict):
            for arg, child in children.iteritems():
                self._children[arg] = child

        # Else, assume it is a list
        else:
            if len(all_arguments) != len(children):
                print "Error: Not all function arguments are mapped to nodes:"
                print "Arguments: ", all_arguments, " Child Nodes: ", children
                raise Exception("Children")
            for arg, child in zip(all_arguments, children):
                self._children[arg] = child

        # Check if there are any undefined nodes for the function
        for arg in all_arguments:
            if arg not in self._children:
                print "Error: Function argument: ", arg,
                print " not given as a node or as a default"
                print "Nodes: ", self._children
                print "Required Function Args: ", all_arguments
                raise Exception()
        pass


    def _evaluate(self):
        kwargs = {}
        for name, node in self._children.iteritems():
            kwargs[name] = node.getVal()
        return self._func(**kwargs)


    def getVal(self):
        if not self._requires_update():
            return self._cached_val
        value = self._evaluate()
        # cache the value
        self._cached_val = value
        for name, child in self._children.iteritems():
            if child.__class__.__name__ == 'variable':
                self._var_values[child.name] = child.val
        self._dirty=False
        return value


    def set_state(self, **kwargs):
        """ Set the state based on the values of the arguments
        Return any args that aren't parameters of the likelihood
        """
        for (arg, val) in kwargs.iteritems():
            self.var(arg).val = val


    def __call__(self, **kwargs):
        # Set values based on kwargs
        # and then return the val
        self.set_state(**kwargs)
        return self.getVal()
    

    def _requires_update(self):
        if self._var_values == {}: 
            self._dirty=True
            return True
        for name, child in self._children.iteritems():
            if child.__class__.__name__ == 'variable': 
                if abs(child.val - self._var_values[child.name]) > self._tolerence:
                    self._dirty=True
                    return True
            else:
                if child._requires_update():
                    self._dirty = True
                    return True
        if self._dirty==True:
            return True

        
    def dependent_vars(self):
        dep_vars = []
        for name, node in self._children.iteritems():
            if node.__class__.__name__ == "variable":
                dep_vars.append(node)
            else:
                dep_vars.extend(node.dependent_vars())
        return dep_vars


    def var(self, var_name):
        for var in self.dependent_vars():
            if var.name == var_name:
                return var
        print "Error: Didn't find variable: ", var_name,
        print " in node: ", self.name
        raise Exception()


    def _numeric_integral(self, vars_to_integrate):
        """ Numerically integrate this node over the given vars

        This will not change the value of the vars after the
        integral has been done (their initial values are restored).

        We specifically handle the case where this node doesn't
        depened on any specified variables by multiplying
        over those variable's ranges.

        """
        
        dep_vars = self.dependent_vars()
        non_dep_vars = [var for var in vars_to_integrate
                        if var not in dep_vars]

        # Handle the case where we integrate over some
        # variables that this node doesn't depend on
        if len(non_dep_vars) != 0:
            dep_vars_to_integrate = [var for var in vars_to_integrate
                                     if var in dep_vars]
            non_dep_integral = 1.0
            for var in non_dep_vars:
                non_dep_integral *= (var.max - var.min)
            if len(dep_vars_to_integrate)==0:
                return non_dep_integral * self.getVal()
            else:
                return non_dep_integral * self._numeric_integral(dep_vars_to_integrate)

        # We have to normalize, but we should be sure
        # to restore the state after, so we aren't effected
        # by the random data points used to evaluate the integral
        vars_before = {var.name : var.val for var in vars_to_integrate}

        # Do the numeric normalization
        if len(vars_to_integrate)==0:
            print "Error: No variables supplied for integration"
            raise Exception("IntegralDimension")

        elif len(vars_to_integrate)==1:
            var_to_integrate = vars_to_integrate[0]
            def func_for_int(var_val):
                var_to_integrate.val = var_val
                return self._evaluate()
            var_min, var_max = (var_to_integrate.min, var_to_integrate.max)
            integral, err = scipy.integrate.quad(func_for_int, var_min, var_max) 
            return integral

        elif len(vars_to_integrate)==2:
            var0 = vars_to_integrate[0]
            var1 = vars_to_integrate[1]
            def func_for_int(var1_val, var0_val):
                var0.val = var0_val
                var1.val = var1_val
                return self._evaluate()
            var0_min, var0_max = (var0.min, var0.max)
            var1_min, var1_max = (var1.min, var1.max)
            integral, err = scipy.integrate.dblquad(func_for_int, var0_min, var0_max,
                                                    lambda x: var1_min, lambda x: var1_max)
            return integral

        elif len(vars_to_integrate)==3:
            var0 = vars_to_integrate[0]
            var1 = vars_to_integrate[1]
            var2 = vars_to_integrate[2]
            def func_for_int(var2_val, var1_val, var0_val):
                var0.val = var0_val
                var1.val = var1_val
                var2.val = var2_val
                return self._evaluate()
            var0_min, var0_max = (var0.min, var0.max)
            var1_min, var1_max = (var1.min, var1.max)
            var2_min, var2_max = (var2.min, var2.max)
            def lower_bound_2(x, y):
                return var2_min
            def upper_bound_2(x, y):
                return var2_max
            integral, err = scipy.integrate.tplquad(func_for_int, var0_min, var0_max,
                                                    lambda x: var1_min, lambda x: var1_max,
                                                    lower_bound_2, upper_bound_2)
            return integral

        
        else: 
            print "Error: Cannot yet integrate over >2 dimensions"
            print "Attempting to integrate over: "
            print [(var, var.name) for var in vars_to_integrate]
            raise Exception("IntegralDimension")


    def integral(self, *vars_to_integrate):
        """ A smart integral that attempts to divide and conquor

        This integration method 
        """

        return self._numeric_integral(vars_to_integrate)


    def __add__(self, other):
        """ Return the sum of two nodes

        """
        name = self.name + "_plus_" + other.name
        return sum_node(name, self, other)


    def __mul__(self, other):
        """ Return the product of two nodes

        """
        name = self.name + "_times_" + other.name
        return product_node(name, self, other)


    def __sub__(self, other):
        """ Return the difference of two nodes

        """
        name = self.name + "_minus_" + other.name
        return diff_node(name, self, other)


class sum_node(node):
    """ A node representing the sum of two nodes

    """

    def __init__(self, name, nodeA, nodeB):
        """ Create a new sum node

        The children are the two supplied nodes,
        and the function is a simple sum of the
        called values of the children nodes.
        """

        def node_sum(a, b):
            return a + b
        node.__init__(self, name, node_sum, [nodeA, nodeB])


class product_node(node):
    """ A node representing the sum of two nodes

    """

    def __init__(self, name, nodeA, nodeB):
        """ Create a new product node

        The children are the two supplied nodes,
        and the function is a simple sum of the
        called values of the children nodes.
        """

        def node_product(a, b):
            return a * b
        node.__init__(self, name, node_product, [nodeA, nodeB])


        def integral(self, *vars_to_integrate):
            """ Return a smart integrate of a product

            There are two cases:
              - If the variables are 

            """

            dep_vars0 = self.children[0].dependent_vars()
            dep_vars1 = self.children[1].dependent_vars()
            
            vars_to_integrate_0_only = [var for var in vars_to_integrate
                                        if var in dep_vars0 
                                        and var not in dep_vars1]
            
            vars_to_integrate_1_only = [var for var in vars_to_integrate
                                        if var in dep_vars1 
                                        and var not in dep_vars0]
            
            vars_for_both = [var for var in vars_to_integrate
                             if var in dep_vars0
                             and var in dep_vars1]

            if len(vars_for_both) != 0:
                return self._numeric_integral(vars_to_integrate)

            else:
                int0 = children[0].integrate(*vars_to_integrate_0_only)
                int1 = children[1].integrate(*vars_to_integrate_1_oly)
                return int0*int1


class diff_node(node):
    """ A node representing the difference of two nodes

    """

    def __init__(self, name, nodeA, nodeB):
        """ Create a new sum node

        The children are the two supplied nodes,
        and the function is a simple sum of the
        called values of the children nodes.
        """

        def node_diff(a, b):
            return a - b
        node.__init__(self, name, node_diff, [nodeA, nodeB])


class quotient_node(node):
    """ A node representing the quotient of two nodes

    """

    def __init__(self, name, nodeA, nodeB):
        """ Create a new quotient node

        The children are the two supplied nodes,
        and the function is a simple sum of the
        called values of the children nodes.
        """

        def node_quotient(a, b):
            return a / b
        node.__init__(self, name, node_quotient, [nodeA, nodeB])



def make_variables(var_string):
    """ Create a tuple of nodes
    from a comma separated list of strings
    """

    # Use a regex 
    # look for: 

    var_list = []

    variables = re.findall(r'(\w+)(\[.*?\])?\s?', var_string)
    for (var, args) in variables:
        #args = re.findall(r'[-]?(\d.)[,]?[\s]?', args)
        #args = re.findall(r'[-]?(\d.)', args)
        #args = re.findall(r'([-]?\d+[.]\d.)', args)
        #args = re.findall(r'([-]?\d+[.]?\d*)', args)

        args = re.findall(r'([.]?[-\d]+[.]?\d*)', args)
        args = map(float, args)
        print "Args: ", args
        if len(args) not in [0,1,3]:
            print "Error: bad arguments for variable %s" % var, args
        var = variable(var, *args)
        var_list.append(var)
    
    return tuple(var_list)

'''
    var = variable(var_name, var_args[0], var_args[1], var_args[2])





    # "x0[.2,-5,5], mu0[0,-5,5], sigma0[1,0,3], fish"

    # Get any variables with ranges defined
    remaining = []
    for var_str in var_string.split('],'):
        if '[' not in var_str and ']' not in var_str:
            remaining.append(var_str)
            print "Appending to remaining: ", var_str
            continue
        var_str = (var_str + ']').strip()
        print "Making ranged variable: ", var_str
        if '[' not in var_str or ']' not in var_str:
            print "Improper variable defined: ", var_str
            raise Exception()
        arg_begin = var_str.find('[')
        arg_end = var_str.find(']')
        if arg_end < arg_begin:
            print "Error, invaid variable string: ", var_str
            raise Exception()
        var_name = var_str[:arg_begin]
        var_args = var_str[arg_begin + 1 : arg_end]
        var_args = var_args.split(',')
        var_args = [ float(arg) for arg in var_args]
        print "Variable args: ", var_args
        if len(var_args)==0:
            print "Error: Improper variable args"
            raise Exception()
        elif len(var_args)==1:
            var = variable(var_name, var_args[0])
        elif len(var_args)==2:
            print "Error: Improper variable args"
            raise Exception()
        elif len(var_args)==3:
            var = variable(var_name, var_args[0], var_args[1], var_args[2])
        else:
            print "Error: Improper variable args"
            raise Exception()
        var_list.append(var)

    for var_str in remaining:
        var_str = var_str.strip()
        var = variable(var_str)
        var_list.append(var)

    return tuple(var_list)
        

'''
