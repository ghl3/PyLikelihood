

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

    '''
    @property
    def val(self):
        return self._val
    @val.setter
    def val(self, val):
        self._val = val
    '''

    def linspace(self):
        return np.linspace(self.min, self.max, self.num_points)

    pass

