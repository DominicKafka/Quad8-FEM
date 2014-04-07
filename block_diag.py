# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 16:13:34 2014

@author: Dominic
"""
# http://comments.gmane.org/gmane.comp.python.scientific.user/20664
import numpy as np

__all__ = ['block_diag']


def block_diag(*arrs):
    """Create a new diagonal matrix from the provided arrays.

    Parameters
    ----------
    a, b, c, ... : ndarray
        Input arrays.

    Returns
    -------
    D : ndarray
        Array with a, b, c, ... on the diagonal.

    """
    arrs = [np.asarray(a) for a in arrs]
    shapes = np.array([a.shape for a in arrs])
    out = np.zeros(np.sum(shapes, axis=0))

    r, c = 0, 0
    for i, (rr, cc) in enumerate(shapes):
        out[r:r + rr, c:c + cc] = arrs[i]
        r += rr
        c += cc
    return out


def test_block_diag():
    x = block_diag(np.eye(2), [[1, 2], [3, 4], [5, 6]], [[1, 2, 3]])
    assert np.all(x == [[1, 0, 0, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0, 0, 0],
                        [0, 0, 1, 2, 0, 0, 0],
                        [0, 0, 3, 4, 0, 0, 0],
                        [0, 0, 5, 6, 0, 0, 0],
                        [0, 0, 0, 0, 1, 2, 3]])