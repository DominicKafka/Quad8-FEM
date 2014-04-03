from numpy import abs

import logging
from scipy.io import loadmat


def build(checkfile):
    checkvalues = loadmat(checkfile, squeeze_me=True)

    def check(name, value, epsilon=1e-10):
        oldvalue = checkvalues[name]
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