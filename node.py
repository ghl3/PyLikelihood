
import inspect
from variable import variable


class node(object):
    """ An node in an acyclic graph

    Can be connected to other nodes to form a graph/tree.
    It uses this tree structure to determine its caching
    and its integration

    """

    def __init__(self, name, func, child_map):

        self.name = name
        self._func = func

        # graph properties
        self._children = {}
        
        # caching properties
        self._var_values = {}
        self._cached_val=None
        self._tolerence = 1E-8
        self._dirty = True

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
        print "Setting Children: ", self._children


    def _evaluate(self):
        kwargs = {}
        for name, node in self._children.iteritems():
            print "Getting value of node: ", name, node, node.__class__.__name__, node.getVal()
            kwargs[name] = node.getVal()
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

    
    def _requires_update(self):
        for name, child in self._children.iteritems():
            if child.__class__.__name__ == 'variable': 
                if self._var_values == {}: 
                    self._dirty=True
                    return True
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

            
