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
        [Fvec[3], 0, 0.5 * Fvec[1], 0.5 * Fvec[1]]])  # Eq.(4.86)
    return Fmat


def Svec_to_Smat(Svec=None):
    import numpy as np
    Smat = np.matrix([[Svec[0], 0, Svec[2], 0], [0, Svec[1], 0, Svec[2]],
        [Svec[2], 0, Svec[1], 0], [0, Svec[2], 0, Svec[0]]])  # Eq.(4.87)
    return Smat