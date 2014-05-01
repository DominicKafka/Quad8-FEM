from numpy import abs

import logging
from scipy.io import loadmat


def build(checkfile):
    checkvalues = loadmat(checkfile, squeeze_me=True)

    def check(name, value, epsilon=5e-8):
        oldvalue = checkvalues[name]
        if (abs(oldvalue - value) > epsilon).any():
            print "'{}' differs from the stored version!".format(name)
            print "Stored value:"
            print oldvalue
            print "Calculated value"
            print value
            print "Difference"
            print oldvalue - value
            raise ValueError
        else:
            (logging.
            debug('{} checked against storedvalue successfully'.format(name)))

    return check