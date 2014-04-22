def B_Quad8(xi, eta, X, NL_flag):
    import numpy as np
    D = np.matrix(np.zeros(4, 16))  # Initialize D with zeros
    # Derivatives of shape functions wrt xi & eta

    dNdxi = np.array([[1 / 4. * eta + 1 / 2. * xi * (1 - eta)
    - 1 / 4 * eta ** 2,
    -1 / 4. * eta + 1 / 2. * xi * (1 - eta) + 1 / 4. * eta ** 2,
    1 / 4. * eta + 1 / 4. * eta ** 2 + 1 / 2. * xi * (1 + eta),
    -1 / 4. * eta + 1 / 2. * xi * (1 + eta) - 1 / 4. * eta ** 2,
    -xi * (1 - eta), 1 / 2. - 1 / 2. * eta ** 2, -xi * (1 + eta),
    -1 / 2. + 1 / 2. * eta ** 2]
    [1 / 4. * xi - 1 / 4. * xi ** 2 + (1 / 2. - 1 / 2. * xi) * eta,
    -1 / 4. * xi - 1 / 4. * xi ** 2 + (1 / .2 + 1 / 2. * xi) * eta,
    1 / 4. * xi + (1 / 2. + 1 / 2. * xi) * eta + 1 / 4. * xi ** 2,
    -1 / 4. * xi + 1 / 4. * xi ** 2 + (1 / 2. - 1 / 2. * xi) * eta,
    -1 / 2. + 1 / 2. * xi ** 2, -2 * (1 / 2. + 1 / 2. * xi) * eta,
    1 / 2. - 1 / 2. * xi ** 2, -2 * (1 / 2. - 1 / 2. * xi) * eta]])
    D[0:1, 0:15:2] = dNdxi
    D[2:3, 2:16:2] = dNdxi  # Arrange shape function derivatives into D
    J = dNdxi * X  # Eq.(2.40)
    detJ = np.linalg.det(J)  # Determinant of Jacobian J
    invJ = np.linalg.inv(J)  # Inverse of Jacobian
    if NL_flag:    # Four rows required for nonlinear case
        dxidX = np.matrix([invJ[0, 0], invJ[0, 1], 0, 0], [0, 0, invJ[1, 0],
             invJ[1, 1]], [invJ[1, 0], invJ[1, 1], 0, 0], [0, 0, invJ[0, 0],
             invJ[0, 1]])
    else:    # Three rows required in linear case
        dxidX = np.matrix([invJ[0, 0], invJ[0, 1], 0, 0], [0, 0, invJ[1, 0],
        invJ[1, 1]], [invJ[1, 0], invJ[1, 1], invJ[0, 0], invJ[0, 1]])

    B = dxidX * D  # Shape function derivatives wrt x and y: Eq.(2.39)
    return B, detJ