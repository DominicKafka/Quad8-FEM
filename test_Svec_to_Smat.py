import unittest
import numpy as np
from solver_utils import Svec_to_Smat


class Svec_to_SmatTestCase(unittest.TestCase):
    def testSvec_to_Smat(self):

        Svec = np.array([[1], [2], [3], [4]])
        Smat = Svec_to_Smat(Svec)

        Smatdash = np.matrix([[1, 0, 3, 0], [0, 2, 0, 3], [3, 0, 2, 0],
        [0, 3, 0, 1]])
        np.testing.assert_almost_equal(Smat,
             Smatdash, decimal=3, err_msg='Incorrect Smat', verbose=True)
        SmatT = np.matrix(Smat).T
        np.testing.assert_almost_equal(SmatT, Smat, 4, 'Not symmetrical')

if __name__ == '__main__':
    unittest.main()