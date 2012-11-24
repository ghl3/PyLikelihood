

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
        return self._val
    @val.setter
    def val(self, val):
        self._val = val

    def linspace(self):
        return np.linspace(self.min, self.max, self.num_points)

    pass


def make_variables(string):
    """ Create a tuple of nodes
    from a comma separated list of strings
    """

    var_list = []
    for var_str in string.split(','):
        var_str = var_name.strip()
        var = None
        if '[' in var_str and ']' in var_str:
            arg_begin = var_str.find(']')
            arg_end = var_str.find('['):
            if arg_end < arg_begin
                print "Error, invaid variable string: ", var_str
                raise Exception()
            var_args = var_str[var_begin : var_end]
            var_args = var_args.split(',')
            if len(var_args)==0:
                print "Error: Improper variable args"
                raise Exception()
            elif len(var_args)==1:
                var = variable(var_name, name=var_args[0])
            elif len(var_args)==2:
                print "Error: Improper variable args"
                raise Exception()
            elif len(var_args)==3:
                var = variable(var_name, var_args[0], var_args[1], var_args[2])
            else:
                print "Error: Improper variable args"
                raise Exception()
        else:
            var = variable(var_name)
        var_list.append(var)
    return tuple(var_list)
        
