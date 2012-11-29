
# cause, duh:
from __future__ import division

import inspect
from variable import variable

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
        print "Setting Children: ", self._children


    def _evaluate(self):
        kwargs = {}
        print "Getting Value of Node: ",
        for name, node in self._children.iteritems():
            print name, node.getVal(),
            kwargs[name] = node.getVal()
        print ''
        return self._func(**kwargs)


    def getVal(self):
        if not self._requires_update():
            print "Using cached val for %s:" % self.name, self._cached_val
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
                    print "Requires update since var: %s changed" % name
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


