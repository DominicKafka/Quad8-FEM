# Basic test from random operation loop in Beam2by20

import numpy as np
import unittest
from scipy.io import loadmat


class SolverUtilsTests(unittest.TestCase):
    def testB_Quad8(self):
        from solver_utils import B_Quad8

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


    def testFvec_to_Fmat(self):
        from solver_utils import Fvec_to_Fmat

        Fvec = np.array([[1],
                         [1],
                         [0],
                         [0]])
        Fmat = Fvec_to_Fmat(Fvec)

        Fmatdash = np.matrix([[1, 0, 0, 0],
                              [0, 1, 0, 0],
                              [0, 0, 0.5, 0.5],
                              [0, 0, 0.5, 0.5]])
        np.testing.assert_almost_equal(Fmat,
             Fmatdash, decimal=3, err_msg='Incorrect Fmat', verbose=True)


    def testQuad8_Res_and_Tangent(self):
        from solver_utils import Quad8_Res_and_Tangent

        # TODO: rewrite this test so that it doesn't load from files
        content = loadmat('Quad8_Res_and_Tangent.mat')
        XY = content['XY']
        U_el = content['U_el']
        matC = content['matC']
        t = content['t']
        El_resd = content['El_res']
        k_elemd = content['k_elem']
        El_stressd = content['El_stress']
        El_straind = content['El_strain']
        [El_res, k_elem, El_stress, El_strain] = Quad8_Res_and_Tangent(XY,
             U_el, matC, t)

        np.testing.assert_almost_equal(El_res,
             El_resd, decimal=3, err_msg='Incorrect El_res', verbose=True)
        np.testing.assert_almost_equal(k_elem,
             k_elemd, decimal=5, err_msg='Incorrect k_elem', verbose=True)
        np.testing.assert_almost_equal(El_stress,
             El_stressd, decimal=3, err_msg='Incorrect El_stress', verbose=True)
        np.testing.assert_almost_equal(El_strain,
             El_straind, decimal=3, err_msg='Incorrect El_strain', verbose=True)


    def testSvec_to_Smat(self):
        from solver_utils import Svec_to_Smat

        Svec = np.array([[1],
                         [2],
                         [3],
                         [4]])
        Smat = Svec_to_Smat(Svec)

        Smatdash = np.matrix([[1, 0, 3, 0],
                              [0, 2, 0, 3],
                              [3, 0, 2, 0],
                              [0, 3, 0, 1]])
        np.testing.assert_almost_equal(Smat,
             Smatdash, decimal=3, err_msg='Incorrect Smat', verbose=True)
        SmatT = np.matrix(Smat).T
        np.testing.assert_almost_equal(SmatT, Smat, 4, 'Not symmetrical')


class MeshTestCase(unittest.TestCase):
    def testNodeCoordDispl(self):
        from mesh_utils import node_coord_displ

        h = 1.
        l = 10.
        m = 10
        n = 100
#        option = 1.
#        E = 200000.
#        nu = 0.3
#        t = 1.
        load_opt = 0
#        Mag = 0
        R = 0
#        RadD = 0
#        V0 = 100.

        coord, displ, elnode, node, el = node_coord_displ(h, l, m, n, load_opt, R)
        self.assertEqual(node, 3221, "Incorrect number of nodes generated")
        # TODO: Fix this check
        self.assertEqual(len(elnode), node, "Elnode is the wrong length")

if __name__ == '__main__':
    unittest.main()
