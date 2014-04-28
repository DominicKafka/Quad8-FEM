import unittest
import numpy as np
from solver_utils import Fvec_to_Fmat


class Fvec_to_FmatTestCase(unittest.TestCase):
    def testFvec_to_Fmat(self):

        Fvec = np.array([[1], [1], [0], [0]])
        Fmat = Fvec_to_Fmat(Fvec)

        Fmatdash = np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0.5, 0.5],
        [0, 0, 0.5, 0.5]])
        np.testing.assert_almost_equal(Fmat,
             Fmatdash, decimal=3, err_msg='Incorrect Fmat', verbose=True)


if __name__ == '__main__':
    unittest.main()