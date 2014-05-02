# -*- coding: utf-8 -*-
import logging


def readsectionfile(filename):
    """ Read a file containing string section names and rows of numbers

    This returns a dictionary with an entry for each section and parses
    the array entries into floats as a list of rows.
    """
    section = {}
    sectionname = None
    with open(filename) as f:
        for line in f:
            try:
                row = [float(item) for item in line.split()]
                section[sectionname].append(row)
            except ValueError:
                sectionname = line.strip()
                section[sectionname] = []
    return section


def read_input_file(filename):
    """ Read sectioned input file for a FEM problem"""

    #TODO: the proper types should be used rather than just float or int for
    #everything
    import time
    from numpy import array

    tic = time.time()

    section = readsectionfile(filename)

    scalarsections = ['Number_of_nodes',
                      'Number_of_elements',
                      'Plane_stress_or_strain',
                      'Number_of_prescribed_displacements',
                      'Number_of_nodal_loads',
                      'Number_of_load_increments',
                      'Number_of_MPCs',
                      ]

    nnodes, nelem, plane, ndispl, ncload, nloadinc, nMPC = (
        [int(section[s][0][0]) for s in scalarsections])
    # TODO: plane should be of type bool

    arraysections = ['Nodal_coordinates',
                     'Element_connectivity',
                     'Nodal_loads',
                     ]

    ndcoor, elnodes, cload = [array(section[s]) for s in arraysections]
    elnodes = array(elnodes, int)
    # TODO: elnodes should be integer

    # sort the rows of displ - keep it as a list
    displ = sorted(section['Prescribed_displacements'])

    nodes = ndcoor[:, 0].tolist()
    coor = ndcoor[:, 1:]

    [[elas, pois, t]] = section['Material_properties']

    MasterDOF = []
    SlaveDOF = []

    toc = time.time()
    time = toc - tic
    (logging.info
    ('Done reading input file {} in {} seconds'.format(filename, time)))
    return (nnodes, ndcoor, nodes, coor, nelem, plane, elnodes, elas, pois, t,
         ndispl, displ, ncload, cload, nloadinc, MasterDOF, SlaveDOF)


def B_Quad8(xi, eta, X, NL_flag):
    import numpy as np
    D = np.matrix(np.zeros([4, 16]))  # Initialize D with zeros
    # Derivatives of shape functions wrt xi & eta

    dNdxi = np.matrix([[1 / 4. * eta + 1 / 2. * xi * (1 - eta)
    - 1 / 4. * eta ** 2,
    -1 / 4. * eta + 1 / 2. * xi * (1 - eta) + 1 / 4. * eta ** 2,
    1 / 4. * eta + 1 / 4. * eta ** 2 + 1 / 2. * xi * (1 + eta),
    -1 / 4. * eta + 1 / 2. * xi * (1 + eta) - 1 / 4. * eta ** 2,
    -xi * (1 - eta), 1 / 2. - 1 / 2. * eta ** 2, -xi * (1 + eta),
    -1 / 2. + 1 / 2. * eta ** 2],
    [1 / 4. * xi - 1 / 4. * xi ** 2 + (1 / 2. - 1 / 2. * xi) * eta,
    -1 / 4. * xi - 1 / 4. * xi ** 2 + (1 / 2. + 1 / 2. * xi) * eta,
    1 / 4. * xi + (1 / 2. + 1 / 2. * xi) * eta + 1 / 4. * xi ** 2,
    -1 / 4. * xi + 1 / 4. * xi ** 2 + (1 / 2. - 1 / 2. * xi) * eta,
    -1 / 2. + 1 / 2. * xi ** 2, -2 * (1 / 2. + 1 / 2. * xi) * eta,
    1 / 2. - 1 / 2. * xi ** 2, -2 * (1 / 2. - 1 / 2. * xi) * eta]])

    for i in range(2):
        for j in range(8):
            D[i, 2 * j] = dNdxi[i, j]
            D[i + 2, 2 * j + 1] = dNdxi[i, j]
              # Arrange shape function derivatives into D
    J = dNdxi * np.matrix(X)  # Eq.(2.40)
    detJ = np.linalg.det(J)  # Determinant of Jacobian J
    invJ = np.linalg.inv(J)  # Inverse of Jacobian
    if NL_flag:    # Four rows required for nonlinear case
        dxidX = np.matrix([[invJ[0, 0], invJ[0, 1], 0, 0], [0, 0, invJ[1, 0],
             invJ[1, 1]], [invJ[1, 0], invJ[1, 1], 0, 0], [0, 0, invJ[0, 0],
             invJ[0, 1]]])
    else:    # Three rows required in linear case
        dxidX = np.matrix([invJ[0, 0], invJ[0, 1], 0, 0], [0, 0, invJ[1, 0],
        invJ[1, 1]], [invJ[1, 0], invJ[1, 1], invJ[0, 0], invJ[0, 1]])

    B = dxidX * D  # Shape function derivatives wrt x and y: Eq.(2.39)
    return B, detJ


def Fvec_to_Fmat(Fvec):
    import numpy as np
    Fmat = np.matrix([[Fvec[0], 0, 0.5 * Fvec[2], 0.5 * Fvec[2]],
        [0, Fvec[1], 0.5 * Fvec[3], 0.5 * Fvec[3]],
        [0, Fvec[2], 0.5 * Fvec[0], 0.5 * Fvec[0]],
        [Fvec[3], 0, 0.5 * Fvec[1], 0.5 * Fvec[1]]], float)  # Eq.(4.86)
    return Fmat


def Svec_to_Smat(Svec=None):
    import numpy as np
    Smat = np.matrix([[Svec[0], 0, Svec[2], 0], [0, Svec[1], 0, Svec[2]],
        [Svec[2], 0, Svec[1], 0], [0, Svec[2], 0, Svec[0]]], float)  # Eq.(4.87)
    return Smat


def Quad8_Res_and_Tangent(X, U, Cmat, t):
    import numpy as np
    Tangent = np.matrix(np.zeros([16, 16]))  # Initialize tangent with zeros
    Res = np.matrix(np.zeros([16, 1]))  # Initialize residual with zeros
    Gauss_pos = 1. / (3) ** (0.5)  # Gauss point location
    Ivec = np.array([[1], [1], [0], [0]])  # Identity vector Eq.(4.89)
    Stress = np.zeros([12, 1])  # Initialize stress vector with zeros
    Strain = np.zeros([12, 1])  # Initialize strain vector with zeros
    true = 1
    Cmat = np.matrix(Cmat)

    for jGauss in range(2):    # 2 by 2 Gauss integration loops
        eta = (-1) ** (jGauss + 1) * Gauss_pos  # Natural coordinate eta
        for iGauss in range(2):
            xi = ((-1) ** (iGauss + 1)) * Gauss_pos   # Natural coordinate xi
            [B, detJ] = B_Quad8(xi, eta, X, true)     # B for Eq.(4.83)
            Fvec = np.array(Ivec + B * np.matrix(U))  # Eq.(4.88)
            Fmat = Fvec_to_Fmat(Fvec)                 # Eq.(4.86)
            Evec = 0.5 * (Fmat.T * Fvec - Ivec)       # Eq.(4.90)
            Svec = Cmat * Evec                        # Eq.(4.96)
            Smat = Svec_to_Smat(Svec)                 # Eq.(4.87)
            Res = Res + B.T * Fmat * Svec * detJ * t  # Eq.(4.98)
            Tangent = (Tangent +
            B.T * (Smat + Fmat * Cmat * Fmat.T) * B * float(detJ) * float(t))
            # Eq.(4.106)
            detF = Fvec[0] * Fvec[1] - Fvec[2] * Fvec[3]
            FTF = np.matrix([
            [float(Fvec[0]) ** 2, Fvec[2] ** 2, Fvec[0] * Fvec[2],
            Fvec[2] * Fvec[0]],
            [Fvec[3] ** 2, Fvec[1] ** 2, Fvec[3] * Fvec[1], Fvec[3] * Fvec[1]],
            [Fvec[0] * Fvec[3], Fvec[2] * Fvec[1], Fvec[0] * Fvec[1],
            Fvec[2] * Fvec[3]]])
            Cauchy = FTF * Svec / detF
            GP = 2 * (jGauss) + iGauss + 1
            for i in range(3):
                f = 3 * (GP - 1) + i
                Stress[f, 0] = Cauchy[i]
                Strain[f, 0] = Evec[i]
    return Res, Tangent, Stress, Strain


def nodal_stresses(elnodes, stress):
    import numpy as np

    StressOut = np.c_[elnodes[:, 0],
         stress[:, [0, 1, 2, 3, 4, 5, 9, 10, 11, 6, 7, 8]]]

    factors = np.matrix([(1 + np.sqrt(3) / 2.), -0.5, -0.5,
         (1 - np.sqrt(3) / 2.)]).T

    cols = np.matrix([[2, 5, 11, 8], [5, 2, 8, 11], [8, 5, 11, 2],
         [11, 2, 8, 5]])

    StressNode = np.matrix(np.zeros(np.shape(StressOut)))
    elnodes = np.matrix(elnodes)
    StressNode[:, 0] = elnodes[:, 0]

    for i in range(3):
        for row in cols:
            row = (row + i - 1)
            StressNode[:, row[0, 0]] = (StressOut[:,
                     [row[0, 0], row[0, 1], row[0, 2], row[0, 3]]] * factors)

    return StressOut, StressNode
