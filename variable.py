

import numpy as np

class variable(object):
    """ A class to store a name, value, and range

    """
    # Constructor
    def __init__(self, name, var_min=-np.inf, var_max=np.inf, num_points=1000):
        self.name = name
        self.min = var_min
        self.max = var_max
        self.num_points = num_points
        self._val = None

    def getVal(self):
        return self._val
    @property
    def val(self):
        print "Getting Val: ", self._val
        return self._val
    @val.setter
    def val(self, val):
        print "Setting Val to: ", val
        self._val = val

    def linspace(self):
        return np.linspace(self.min, self.max, self.num_points)

    pass


def make_variables(string):
    """ Create a tuple of nodes
    from a comma separated list of strings
    """

    var_list = []
    for var_name in string.split(','):
        var_name = var_name.strip()
        print "Making Var with name: ", var_name
        var = variable(var_name)
        var_list.append(var)
    return tuple(var_list)
        
