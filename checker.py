from numpy import abs

import logging

def build(checkvalues):

    def check(name, value, epsilon=1e-10):
        oldvalue = checkvalues[name]
        # We need to do some checking of dimensions to compare correctly
        if oldvalue.shape[0] == 1 or oldvalue.shape[1] == 1: # vector
            oldvalue = oldvalue.flatten()
        if (abs(oldvalue - value) > epsilon).any():
            print "'{}' differs from the stored version!".format(name)
            print "Stored value:"
            print oldvalue
            print "Calculated value"
            print value
            raise ValueError
        else:
            logging.debug('{} checked against stored value successfully'.format(name))

    return check