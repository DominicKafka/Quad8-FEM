# -*- coding: utf-8 -*-
import logging

import numpy as np

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
                      'Number_of_height_increments',
                      'Number_of_length_increments',
                      ]

    nnodes, nelem, plane, ndispl, ncload, nloadinc, nMPC, m, n = (
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
    thickness = section['Spring_thickness'][0][0]

    MasterDOF = []
    SlaveDOF = []

    # Check for consistency
    assert nnodes == len(ndcoor)
    assert nelem == len(elnodes)
    assert ndispl == len(displ)
    assert ncload == len(cload)

    toc = time.time()
    time = toc - tic
    (logging.info
    ('Done reading input file {} in {} seconds'.format(filename, time)))
    return (ndcoor, nodes, coor, plane, elnodes, elas, pois, t, displ, cload,
         nloadinc, m, n, thickness, MasterDOF, SlaveDOF)


def B_Quad8(xi, eta, X, NL_flag):

    D = np.asmatrix(np.zeros([4, 16]))  # Initialize D with zeros
    # Derivatives of shape functions wrt xi & eta

    dNdxi = 1/4.*np.array([[eta + 2. * xi * (1 - eta) - eta ** 2,
                             -eta + 2. * xi * (1 - eta) + eta ** 2,
                             eta + eta ** 2 + 2. * xi * (1 + eta),
                             -eta + 2. * xi * (1 + eta) - eta ** 2,
                             4*(-xi * (1 - eta)),
                             2. - 2. * eta ** 2,
                             4*(-xi * (1 + eta)),
                             -2. + 2. * eta ** 2],
                            [xi - xi ** 2 + (2. - 2. * xi) * eta,
                             -xi - xi ** 2 + (2. + 2. * xi) * eta,
                             xi + (2. + 2. * xi) * eta + xi ** 2,
                             -xi + xi ** 2 + (2. - 2. * xi) * eta,
                             -2. + 2. * xi ** 2,
                             -4 * (1 + xi) * eta,
                             2. - 2. * xi ** 2,
                             -4 * (1 - xi) * eta]])

    # Arrange shape function derivatives into D
    i = np.arange(2)
    j = np.arange(8)
    D[np.ix_(i, 2 * j)] = dNdxi
    D[np.ix_(i + 2, 2 * j + 1)] = dNdxi

    J = dNdxi * np.matrix(X)  # Eq.(2.40)
    detJ = np.linalg.det(J)  # Determinant of Jacobian J
    invJ = np.linalg.inv(J)  # Inverse of Jacobian
    if NL_flag:    # Four rows required for nonlinear case
        dxidX = np.matrix([[invJ[0, 0], invJ[0, 1], 0, 0],
                           [0, 0, invJ[1, 0], invJ[1, 1]],
                           [invJ[1, 0], invJ[1, 1], 0, 0],
                           [0, 0, invJ[0, 0], invJ[0, 1]]])
    else:    # Three rows required in linear case
        dxidX = np.matrix([invJ[0, 0], invJ[0, 1], 0, 0],
                          [0, 0, invJ[1, 0], invJ[1, 1]],
                          [invJ[1, 0], invJ[1, 1], invJ[0, 0], invJ[0, 1]])

    B = dxidX * D  # Shape function derivatives wrt x and y: Eq.(2.39)
    return B, detJ


def Fvec_to_Fmat(Fvec):

    Fmat = np.matrix([[Fvec[0], 0, 0.5 * Fvec[2], 0.5 * Fvec[2]],
                      [0, Fvec[1], 0.5 * Fvec[3], 0.5 * Fvec[3]],
                      [0, Fvec[2], 0.5 * Fvec[0], 0.5 * Fvec[0]],
                      [Fvec[3], 0, 0.5 * Fvec[1], 0.5 * Fvec[1]]], float)  # Eq.(4.86)
    return Fmat


def Svec_to_Smat(Svec=None):

    Smat = np.matrix([[Svec[0], 0, Svec[2], 0],
                      [0, Svec[1], 0, Svec[2]],
                      [Svec[2], 0, Svec[1], 0],
                      [0, Svec[2], 0, Svec[0]]], float)  # Eq.(4.87)
    return Smat


def Quad8_Res_and_Tangent(X, U, Cmat, t):

    Tangent = np.asmatrix(np.zeros([16, 16]))  # Initialize tangent with zeros
    Res = np.asmatrix(np.zeros([16, 1]))  # Initialize residual with zeros
    Gauss_pos = 1. / (3) ** (0.5)  # Gauss point location
    Ivec = np.array([[1], [1], [0], [0]])  # Identity vector Eq.(4.89)
    Stress = np.zeros([12, 1])  # Initialize stress vector with zeros
    Strain = np.zeros([12, 1])  # Initialize strain vector with zeros
    Cmat = np.asmatrix(Cmat)
    i = np.arange(3)
    FTF = np.empty((3, 4))
    for jGauss in range(2):    # 2 by 2 Gauss integration loops
        eta = (-1) ** (jGauss + 1) * Gauss_pos  # Natural coordinate eta
        for iGauss in range(2):
            xi = ((-1) ** (iGauss + 1)) * Gauss_pos   # Natural coordinate xi
            [B, detJ] = B_Quad8(xi, eta, X, True)     # B for Eq.(4.83)
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
            F = Fvec  # PEP8 note: This is more legible than the longer lines
            # Note: Updating FTF in-place like this saves time
            FTF[:] = [[float(F[0])**2, F[2]**2, F[0]*F[2], F[2]*F[0]],
                      [F[3]**2, F[1]**2, F[3]*F[1], F[3]*F[1]],
                      [F[0]*F[3], F[2]*F[1], F[0]*F[1], F[2]*F[3]]]
            Cauchy = FTF * Svec / detF
            GP = 2 * (jGauss) + iGauss + 1
            f = 3 * (GP - 1) + i
            Stress[f, :] = Cauchy[i, 0]
            Strain[f, :] = Evec[i, 0]
    return Res, Tangent, Stress, Strain


def nodal_stresses(elnodes, stress):
    StressOut = np.c_[elnodes[:, 0],
         stress[:, [0, 1, 2, 3, 4, 5, 9, 10, 11, 6, 7, 8]]]

    factors = np.matrix([(1 + np.sqrt(3) / 2.), -0.5, -0.5,
         (1 - np.sqrt(3) / 2.)]).T

    cols = np.matrix([[2, 5, 11, 8],
                      [5, 2, 8, 11],
                      [8, 5, 11, 2],
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


def calc_von_mises(StressNode, pois, plane):
    VonMises = np.array(np.zeros([len(StressNode), 4]))
    StressNode = np.array(StressNode)
    if (plane == 1):
        for i in range(4):
            VonMises[:, i] = (StressNode[:, 1 + i * 3] ** 2 -
                StressNode[:, 1 + i * 3] * StressNode[:, 2 + i * 3] +
                StressNode[:, 2 + i * 3] ** 2 +
                3 * StressNode[:, 3 + i * 3] ** 2)
            VonMises[:, i] = VonMises[:, i] ** 0.5

    else:
        for i in range(4):
            VonMises[:, i] = ((1 - pois + pois ** 2) *
                (StressNode[:, 1 + i * 3] ** 2 + StressNode[:, 2 + i * 3] ** 2)
                 - (1 + pois - pois ** 2) * StressNode[:, 1 + i * 3] *
                 StressNode[:, 2 + i * 3] + 3 * StressNode[:, 3 + i * 3] ** 2)
            VonMises[:, i] = VonMises[:, i] ** 0.5

    return VonMises


def calc_tresca(StressNode, pois, plane):
    nelem = len(StressNode)

    Tresca = np.matrix(np.zeros([nelem, 4]))

    if plane:
        for j in range(nelem):
            for i in range(4):
                s = [[StressNode[j, 1 + i * 3], StressNode[j, 3 + i * 3], 0],
                     [StressNode[j, 3 + i * 3], StressNode[j, 2 + i * 3], 0],
                     [0, 0, 0]]
                principal = np.linalg.eig(s)
                Tresca[j, i] = np.max(principal[0]) - np.min(principal[0])
    else:
        for j in range(nelem):
            for i in range(4):
                s = [[StressNode[j, 1 + i * 3], StressNode[j, 3 + i * 3], 0],
                     [StressNode[j, 3 + i * 3], StressNode[j, 2 + i * 3], 0],
                     [0, 0, pois * (StressNode[j, 1 + i * 3] + StressNode[j, 2 + i * 3])]]
                principal = np.linalg.eig(s)
                Tresca[j, i] = np.max(principal[0]) - np.min(principal[0])

    return Tresca


def write_section(fid, title, data, headings, formats):
    """ Write a section of output file using the title supplied and matching centered
        headings and formats """
    def writeline(string):
        fid.write(string + ' \n')

    # Figure out the total width assuming formats look like %0.0f or %0f
    widths = [int(f[1:-1].split('.')[0]) for f in formats]
    totalwidth = sum(widths) + len(formats) - 1

    writeline(title.center(min(totalwidth, 42)))
    stars = '*' * totalwidth
    writeline(stars)
    writeline(' '.join(h.center(w) for h, w in zip(headings, widths)))
    writeline(stars)
    for row in data:
        writeline(' '.join(formats) % tuple(row))
    writeline('')


def write_output_file(file_out, U, displ, Pb, nodes, elnodes, strain,
    StressOut):

    Uoutput = np.c_[nodes, U[::2], U[1::2]]
    Uoutput = np.array(Uoutput)
    StressOut = np.array(StressOut)

    StrainOut = np.c_[elnodes[:, 0],
        strain[:, [0, 1, 2, 3, 4, 5, 9, 10, 11, 6, 7, 8]]]
    StrainOut = np.array(StrainOut)

    displ = np.array(displ)
    SupReac = np.c_[displ[:, [0, 1]], Pb]
    SupReac = np.array(SupReac)

    fid = open(file_out, 'w')

    fid.write('OUTPUT OF PYTHON Q4 SMALL STRAIN FEM IMPLEMENTATION \n')
    fid.write('\n')
    write_section(fid, 'DISPLACEMENTS', Uoutput,
                  ['Node', 'U1', 'U2'],
                  ['%5d', '%13.5e', '%13.5e'])

    write_section(fid, 'ELEMENT STRESSES', StressOut,
                  ['Element', 'S11_G1', 'S22_G1', 'S12_G1',
                              'S11_G2', 'S22_G2', 'S12_G2',
                              'S11_G3', 'S22_G3', 'S12_G3',
                              'S11_G4', 'S22_G4', 'S12_G4'],
                  ['%5d'] + ['%12.4e']*12)

    write_section(fid, 'ELEMENT STRAINS', StrainOut,
                  ['Element', 'E11_G1', 'E22_G1', 'E12_G1',
                              'E11_G2', 'E22_G2', 'E12_G2',
                              'E11_G3', 'E22_G3', 'E12_G3',
                              'E11_G4', 'E22_G4', 'E12_G4'],
                  ['%5d'] + ['%12.4e']*12)

    write_section(fid, 'SUPPORT REACTIONS', SupReac,
                  ['Node', 'DOF', 'Magnitude'],
                  ['%5d', '%5d', '%17.5e'])
    fid.close()


def graphs(nnodes, coor, nelem, elnodes, StressNode, U,
    VonMises, Tresca, nloadinc, All_soln):
    import numpy as np
    import pylab as pl
    inp = '9'
    while inp.strip() != '0':
        inp = raw_input("Please choose a graph to display:\n 1: Original Shape\n 2: Deformed Shape\n 3: Overlay\n 0: End\n")
        if inp.strip() == '1':
            pl.plot(coor[:, 0], coor[:, 1], 'b.')
            pl.axis('equal')
            pl.show()
        elif inp.strip() == '2':
            dpm = np.c_[U[0:2 * nnodes:2], U[1:2 * nnodes:2]]
            magfac = 1.
            dcoor = coor + magfac * dpm
            pl.plot(dcoor[:, 0], dcoor[:, 1], 'r.')
            pl.axis('equal')
            pl.show()
        elif inp.strip() == '3':
            dpm = np.c_[U[0:2 * nnodes:2], U[1:2 * nnodes:2]]
            magfac = 1.
            dcoor = coor + magfac * dpm
            pl.plot(coor[:, 0], coor[:, 1], 'b.',
                    dcoor[:, 0], dcoor[:, 1], 'r.')
            pl.axis('equal')
            pl.show()


def linear_estimation():

    import pylab as pl
    import numpy as np

    l = 0.64
    b = 0.045
    h = 0.006
    E = 210000000000.
    d1 = 0.01
    d2 = 0.02
    d3 = 0.03
    factor1 = 1.6
    factor2 = 2.2


    I = float(1./12*b*h**3)
    knat = ((48.0*E*I)/(l**3))
    smodulus = E * I

    k1 = knat*factor1
    k2 = knat*factor2

    # Natural Support
    disp1 = np.linspace(0, d1, 10)
    Fx1 = knat * disp1
    F1 = d1*knat
    stress1 = ((0.5*F1*0.5*l*0.5*h)/I)/1000000


    # First Point
    l1 = ((48.*E*I)/k1)**(1./3.)
    x1 = 0.5*l-0.5*l1
    y1 = (F1*x1)/(48.*E*I)*(4.*x1**2-3.*l**2)

    Fcoeff = (((l1 **2) * x1)/(16 * E * I))
    Facoeff = (((x1 ** 2) * l1)/(2. * E * I))+((x1 ** 3) / (3. * E * I))
    Fa1 = ((-y1))/(Facoeff)
    Fcoeff2 = Fcoeff/Facoeff
    Fcoeff3 = (l1**3)/(48*E*I)
    gencoeff = (x1*(l1**2))/(8.*E*I)
    Fcoeff4 = Fcoeff2*gencoeff
    const = Fa1 * gencoeff
    Fcoeff5 = Fcoeff3-Fcoeff4
    Fconst = const/Fcoeff5

    delta1 = (y1-(Fconst*Fcoeff5)-(y1*k1*Fcoeff5))/((k1*Fcoeff5)-1)
    disp12 = np.linspace(d1, delta1, 10)

    Fx121 = (disp12+y1)/Fcoeff5 -Fconst

    Fa = Fa1-((k1*(delta1+y1))*Fcoeff2)
    Faint1 = Fa1-((knat*(0.01))*Fcoeff2)

    # Line 2
    disp2 = np.linspace(delta1, d2, 10)
    #disp2 = np.linspace(0, d2, 10)
    Fx2 = k1 * (disp2-(-y1))
    F2 = k1 * (d2-(-y1))
    stress2 = ((0.5*F2*0.5*l1*0.5*h)/I)/1000000


    # Second Point
    l2 = ((48.*E*I)/k2)**(1./3.)
    x2 = 0.5*l-0.5*l2
    y2 = y1 + (F2*(x2-x1))/(48.*E*I)*(4.*(x2-x1)**2-3.*l1**2)
    ymax2 = y1 + (-1*(F2*l2**3.)/(48.*E*I))

    x12 = x2-x1
    F2coeff = (((l2 **2) * x12)/(16 * E * I))
    Facoeff2 = (((x12 ** 2) * l2)/(2. * E * I))+((x12 ** 3) / (3. * E * I))
    Fa12 = ((-y2-(-y1)))/(Facoeff2)
    Fcoeff22 = F2coeff/Facoeff2
    Fcoeff32 = (l2**3)/(48*E*I)
    gencoeff2 = (x12*(l2**2))/(8.*E*I)
    Fcoeff42 = Fcoeff22*gencoeff2
    const2 = Fa12 * gencoeff2
    Fcoeff52 = Fcoeff32-Fcoeff42
    Fconst2 = const2/Fcoeff52

    delta12 = (y2-(Fconst2*Fcoeff52)-(y2*k2*Fcoeff52))/((k2*Fcoeff52)-1)
    disp122 = np.linspace(d2, delta12, 10)

    Fx1212 = (disp122+y2)/Fcoeff52 -Fconst2

    Fa2 = Fa12-((k2*(delta12+y2))*Fcoeff22)
    Faint12 = Fa12-((k1*(0.02+y1))*Fcoeff22)


    disp3 = np.linspace(delta12, d3, 10)
    #disp3 = np.linspace(0, d3, 10)
    Fx3 = k2 * (disp3+y2)
    F3 = k2 * (d3+y2)
    stress3 = ((0.5*F3*0.5*l2*0.5*h)/I)/1000000


    print 'Moment of inertia:'+str(I)
    print 'Stress1: '+str(stress1)+' MPa'
    print 'x1: '+str(x1)+' m'
    print 'y1: '+str(y1)+' m'
    print l2
    print x2
    print 'Stress2: '+str(stress2)+' MPa'
    print 'x2: '+str(x2)+' m'
    print 'y2: '+str(y2)+' m'
    print 'Stress3: '+str(stress3)+' MPa'

    return disp1, Fx1, disp2, Fx2, disp3, Fx3, disp122, Fx1212, disp12, Fx121


def springtest8():
    import csv
    import numpy as np
    import pylab as pl

    displacement = []
    Force = []
    displacement1 = []
    Force1 = []
    absdisp = []
    rowcount = 0
    absoluted = 0
    absoluteF = 0
    absF = []
    LCfactor = 2 * 9.81 * 1000 / 10
    displacement2 = []
    Force2 = []
    displacement3 = []
    Force3 = []

    reader = csv.reader(open("springtest7.csv", "rb"), delimiter=',')
    x = list(reader)
    result = np.array(x)

    [x, y] = np.shape(result)
    print 'number of rows ' + str(x)
    print 'number of columns ' + str(y)

    for i in range(x):
        if i > 0:
            displacement.append(result[i, 1])
            Force.append(result[i, 4])


    maxdisp = max(displacement)
    mindisp = min(displacement)

    maxF = max(Force)
    minF = min(Force)
    # print Force


    for k in range(x - 1):
        absoluted = float(maxdisp) - float(displacement[k])
        absdisp.append(absoluted)
        absoluteF = (float(minF) - float(Force[k])) * LCfactor
        absF.append(absoluteF)

    for j in range(640):
        if j > 60:
            displacement1.append(absdisp[j])
            Force1.append(absF[j])

    print 'maxdisp ' + str(maxdisp)
    print 'mindisp ' + str(mindisp)
    print 'maxF ' + str(maxF)
    print 'minF ' + str(minF)
    displacement1 = [x * 0.001 for x in displacement1]
    displacement1 = [x + 0.0035 for x in displacement1]
    Force1 = [x + 130 for x in Force1]

    return displacement1, Force1