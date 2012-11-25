

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


def make_variables(var_string):
    """ Create a tuple of nodes
    from a comma separated list of strings
    """

    var_list = []

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
        
