
import inspect
from variable import variable


class node(object):
    """ An node in an acyclic graph

    Can be connected to other nodes to form a graph/tree.
    It uses this tree structure to determine its caching
    and its integration

    """

    _name = ""
    _func = None

    # graph properties
    _children = {}
    #_direct_variables = []
    #_parents = []
    #_child_name_map = {}

    # caching properties
    _var_values = {}
    _cached_val=None
    _tolerence = 1E-8
    _dirty = True

    def __init__(self, name, func, child_map):
        self._name = name
        func_spec = inspect.getargspec(func)
        (all_arguments, all_defaults) = (func_spec.args, func_spec.defaults)
        '''
        # Get any function arguments that
        # have default values
        default_dict = {}
        for name, val in zip(reversed(all_arguments), reversed(all_defaults)):
            default_dict[name] = val
        '''
        # Add the child nodes as children to this node
        for child, node in child_map.iteritems():
            self._children[child] = node
        # Check if there are any undefined nodes for the function
        for arg in all_arguments:
            if arg not in self._children:
                print "Error: Function argument: ", arg,
                print " not given as a node or as a default"
                print "Nodes: ", self._children
                print "Required Function Args: ", all_arguments
                raise Exception()
        self._func = func

    def _evaluate(self):
        kwargs = {}
        for name, node in self._children.iteritems():
            kwargs[name] = node.val
        return self._func(**kwargs)

    def val(self):
        if not self._requires_update():
            return self._cached_val
        value = self._evaluate()
        # cache the value
        self._cached_val = value
        for name, child in self._children.iteritems():
            if child.__class__.__name__ == 'variable':
                self._var_values[child.name] = child.val
        return value
    
    def _requires_update(self):
        for name, child in self._children.iteritems():
            print "Child: ", name, child, child.__class__.__name__
            if child.__class__.__name__ == 'variable': 
                if self._var_values == {}: 
                    self._dirty=True
                    return True
                if abs(child.val - self._var_values[child.name]) > self._tolerence:
                    self._dirty=True
                    return True
            if child._dirty:
                self._dirty = True
                return True
                
