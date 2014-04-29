import unittest
import numpy as np
from solver_utils import Quad8_Res_and_Tangent
from scipy.io import loadmat


class Quad8_Res_and_TangentTestCase(unittest.TestCase):
    def testQuad8_Res_and_Tangent(self):

        content = loadmat('Quad8_Res_and_Tangent.mat')
        #print content
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
             k_elemd, decimal=3, err_msg='Incorrect k_elem', verbose=True)
        np.testing.assert_almost_equal(El_stress,
             El_stressd, decimal=3, err_msg='Incorrect El_stress', verbose=True)
        np.testing.assert_almost_equal(El_strain,
             El_straind, decimal=3, err_msg='Incorrect El_strain', verbose=True)

if __name__ == '__main__':
    unittest.main()