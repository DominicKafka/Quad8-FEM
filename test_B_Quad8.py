# Basic test from random operation loop in Beam2by20

import unittest
import numpy as np
from solver_utils import B_Quad8


class B_Quad8TestCase(unittest.TestCase):
    def testB_Quad8(self):

        xi = 0.57735
        eta = -0.57735
        X = np.array([[8, -0.5], [8.5, -0.5], [8.5, 0.], [8, 0], [8.25, -0.5],
        [8.5, -0.25], [8.25, 0.], [8, -0.25]])
        true = 1

        [B, detJ] = B_Quad8(xi, eta, X, true)

        Bdash = np.matrix([[0.91068, 0.0000, 2.73205, 0.0000, 0.24402, 0.0000,
         0.73205, 0.0000, -3.64273, 0.0000, 1.33333, 0.0000, -0.97607, 0.0000,
          -1.33333, 0.0000],
        [0.0000, -0.24402, 0.0000, -2.73205, 0.0000, -0.91068, 0.0000, -0.73205,
         0.0000, -1.33333, 0.0000, 3.64273, 0.0000, 1.33333, 0.0000, 0.97607],
        [-0.24402, 0.0000, -2.73205, 0.0000, -0.91068, 0.0000, -0.73205, 0.0000,
         -1.33333, 0.0000, 3.64273, 0.0000, 1.33333, 0.0000, 0.97607, 0.0000],
    [0.0000, 0.91068, 0.0000, 2.73205, 0.0000, 0.24402, 0.0000, 0.73205, 0.0000,
         -3.64273, 0.0000, 1.33333, 0.0000, -0.97607, 0.0000, -1.33333]])

        detJdash = 0.062500

        np.testing.assert_almost_equal(B,
             Bdash, decimal=3, err_msg='Incorrect B Matrix', verbose=True)

        np.testing.assert_almost_equal(detJ,
             detJdash, 4, "Incorrect determinant")

if __name__ == '__main__':
    unittest.main()