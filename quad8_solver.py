# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 09:39:32 2014

@author: Dominic
"""
import math
import numpy as np
from scipy.sparse import coo_matrix
from solver_utils import read_input_file
from solver_utils import Quad8_Res_and_Tangent
from block_diag import block_diag
import time

import logging
logging.basicConfig(level=logging.DEBUG)

import checker

#ElType = '5B';
#ElType = '7B';
#ElType = 'q4';
ElType = 'q8'

NodesPerEl = 0
if ElType in ['q4', '5B', '7B']:
    NodesPerEl = 4
elif ElType == 'q8':
    NodesPerEl = 8

DofPerEl = 2 * NodesPerEl
TotKvec = DofPerEl ** 2

print 'Welcome to the Quad-8 Finite Element Program'

#filename = raw_input('Input filename (without extension) ? ')
#print filename
case = 'Beam2by20'

check = checker.build(case + '.mat')
#check = checker.build(case + '.inp')

filename = case + '.inp'
([nnodes, ndcoor, nodes, coor, nelem, plane, elnodes, elas, pois, t, ndispl,
displ, ncload, cload, nloadinc, mdof, sdof]) = read_input_file(filename)
file_out = filename + '.out'
sdof = np.array(sdof)
mdof = np.array(mdof)
check('elnodes', elnodes)

#Get maximum dimensions of model
dx_max = max(coor[:, 0]) - min(coor[:, 0])
dy_max = max(coor[:, 1]) - min(coor[:, 1])
dL_max = math.sqrt(dx_max ** 2 + dy_max ** 2)

#GraphOpt=1: graphical display of results and text based output files
#GraphOpt=0: only text based output files
GraphOpt = 0

# Find prescribed (pdof) and free (fdof) degrees of freedom
# FIXME: the next line assumes 2d
dof = np.ones([nnodes, 2])
Up = []
for node, dimension, displacement in displ:
    dof[nodes.index(node), dimension - 1] = 0
    Up.append(displacement)
Up = np.matrix(Up).T
pdof = np.flatnonzero(dof == 0)
fdof = np.flatnonzero(dof != 0)
#TODO: check if fdof and pdof should actually be a set instead of an array
fdof = np.setdiff1d(fdof, mdof)
fdof = np.setdiff1d(fdof, sdof)

check('fdof', fdof + 1)
check('pdof', pdof + 1)

# Initially guess that all free displacements are zero
U = np.zeros([2 * nnodes, 1])

# Construct elasticity tensor C
# If plane strain
if plane == 0:
    e = elas / float(1 - pois ** 2)
    nu = pois / float(1 - pois)
else:
    e = elas
    nu = pois

c = e / float(1 - nu ** 2)
matC = (np.asmatrix(
    block_diag([[c, c * nu], [c * nu, c]], np.eye(2) * c * (1 - nu))))
check('matC', matC)

ndnum = range(1, (1 + NodesPerEl))
[colpos, rowpos] = np.meshgrid(range(DofPerEl), range(DofPerEl))

colpos = colpos.flatten()
rowpos = rowpos.flatten()


tol = 3e-5
dUNrm = 1.0

F = np.zeros([2 * nnodes, 1])
LoadFac = (np.array(range(1, nloadinc + 1))) / float(nloadinc)
check('LoadFac', LoadFac)

# FIXME: this 12 should be a variable
stress = np.matrix(np.zeros([nelem, 12]))
strain = np.matrix(np.zeros([nelem, 12]))

AllResNrm = []
AlldUNrm = []
All_iter = []
All_soln = []

for iter_load in range(nloadinc):
    print '---------------------------------------------------------------'
    print '                      Load increment ' + str(iter_load + 1)
    print '---------------------------------------------------------------'
    # Correct fraction of prescribed displacement applied to pdof
    U[pdof, :] = LoadFac[iter_load] * Up
    # If MPCs are present, compute slave DOF values
    if (sdof.size == True):
        [U[sdof, 0], P_mpc, d2PdUm2_Rs] = MPC_user(U[mdof, 0], F[sdof])
    else:
        P_mpc = np.matrix([])
        d2PdUm2_Rs = []

    ResNrm = 1.

    itera = 0

    while (ResNrm > tol) or (dUNrm > tol):

        tic = time.time()
        itera = itera + 1
        #Main loop over elements. Compute k_elem and assemble
        #Initialize global stiffness matrix vectors;
        row_vec = np.matrix(np.zeros([TotKvec * nelem, 1]))
        col_vec = np.matrix(np.zeros([TotKvec * nelem, 1]))
        stiff_vec = np.matrix(np.zeros([TotKvec * nelem, 1]))
        Residual = np.matrix(np.zeros([2 * nnodes, 1]))
        #Initialize global load vector
        F_ext = np.matrix(np.zeros([2 * nnodes, 1]))
        pos_vec = 0
        for i in range(nelem):
            # Find reference coordinates of element nodes
            X = coor[elnodes[i, ndnum] - 1, 0]
            Y = coor[elnodes[i, ndnum] - 1, 1]
            # Get global degree of freedom numbers per element
            pg = np.matrix(np.zeros([DofPerEl, 1], int))
            pg[::2, 0] = np.matrix(2 * elnodes[i, 1:(1 + NodesPerEl)] - 2).T
            pg[1::2, 0] = np.matrix(2 * elnodes[i, 1:(1 + NodesPerEl)] - 1).T
            # Get current guess for nodal displacements
            U_el = np.matrix(U[pg], float).T
            XY = np.matrix([X, Y]).T


            if ElType == 'q4':
                [El_res, k_elem, El_stress, El_strain] = Quad4_Res_and_Tangent(
                    XY, U_el, matC, t)
            elif ElType == '5B':
                [El_res, k_elem, El_stress, El_strain] = FiveB_Res_and_Tangent(
                    XY, U_el, matC, t)
            elif ElType == '7B':
                [El_res, k_elem, El_stress, El_strain] = (
                    SevenB_Res_and_Tangent(XY, U_el, matC, t))
            elif ElType == 'q8':
                [El_res, k_elem, El_stress, El_strain] = Quad8_Res_and_Tangent(
                    XY, U_el, matC, t)
            stress[i, :] = El_stress.T
            strain[i, :] = El_strain.T

            # Assemble residual
            Residual[pg, 0] = Residual[pg, 0] + El_res

            # Assemble k_elem into sparse k_global using vectors
            k = TotKvec * i + np.matrix(range(TotKvec)).T
            row_vec[k, 0] = pg[rowpos, 0]
            col_vec[k, 0] = pg[colpos, 0]
            k_elem = np.matrix(k_elem.flatten()).T
            stiff_vec[k, 0] = k_elem[0:TotKvec, 0]

            if (i == 4) and (itera == 1) and (iter_load == 0):
                check = checker.build('nelemloop.mat')
                check('shape_row_vec', np.shape(row_vec))
                check('shape_col_vec', np.shape(col_vec))
                check('shape_stiff_vec', np.shape(stiff_vec))
                #check('row_vec', row_vec)  # gives 'memory error'
                #check('col_vec', col_vec)  # gives 'memory error'
                #check('k_elem', k_elem)  # gives 'memory error'
                #check('stiff_vec', stiff_vec) # gives 'memory error'
                #check('some_stiff_vec', some_stiff_vec.T)
                #gives error though identical
                check('U_el', U_el)
                check('XY', XY)
                check('stress', stress)
                check('strain', strain)
                check('Residual', Residual)
            # End of main loop over elements

        # Assemble k_global from vectors
        k_global = coo_matrix((stiff_vec.T.tolist()[0], (row_vec.T.tolist()[0],
             col_vec.T.tolist()[0])), shape=(2 * nnodes, 2 * nnodes))
        k_global = k_global.todense().T
        toc = time.time()
        finish = toc - tic
        print 'Done assembling stiffness matrix:', finish, 'seconds.'

        # Add nodal loads to global load vector
        for i in range(ncload):
            p = np.where(nodes == cload[i, 0])
            pos = (p[0]) * 2 + cload[i, 1] - 1
            pos = np.array(pos, int)
            bob = F_ext[pos, 0]
            ans = F_ext[pos, 0] + LoadFac[iter_load] * cload[i, 2]
            F_ext[pos, 0] = ans[0, 0]

            if (i == 4) and (itera == 1) and (iter_load == 0):
                check = checker.build('ncloadloop.mat')
                check('F_ext', F_ext.T)

        # Subtract internal nodal loads
        F = Residual - F_ext

        ResNrm = np.linalg.norm(F[fdof, 0])
        if itera == 1:
            if ResNrm > 1e-4:
                ResNrm0 = ResNrm
            else:
                ResNrm0 = 1

        ResNrm = ResNrm / ResNrm0
        print ('Normalized residual at start of iteration ', itera, ' = ',
        ResNrm)
        # Solve update of free dof's
        # Solution for non-symmetric stiffness matrix

        Kff = k_global[np.ix_(fdof, fdof)]
        Pf = F[fdof]
        # Define RHS
        if (mdof.size > 0):
            Kfm = (k_global[np.ix_(fdof, mdof)] +
                k_global[np.ix_(fdof, sdof)] * P_mpc)
            Kmm = (k_global[np.ix_(mdof, mdof)] + k_global[np.ix_(mdof, sdof)]
                * P_mpc + P_mpc.T * k_global[np.ix_(sdof, mdof)] +
                 P_mpc.T * k_global[np.ix_(sdof, sdof)] * P_mpc + d2PdUm2_Rs)
            Pm = F[mdof] + P_mpc.T * (F[sdof])
            Pa = [Pf, Pm]
            Kaa = np.matrix([Kff, Kfm], [Kfm.T, Kmm])
        else:
            Kaa = Kff
            Pa = Pf

        toc = time.time()
        finish = toc - tic
        print 'Done assembling stiffness matrix:', finish, 'seconds.'
        tic = time.time()

        deltaUf = -1 * np.linalg.solve(Kaa, Pa)

        finish = time.time() - tic
        print 'Done solving system:', finish, 'seconds.'

        dUNrm = np.linalg.norm(deltaUf) / dL_max
        print 'Normalized displacement update                 = ', dUNrm
        print '                    --------------------'

        # Sort Uf and Ub into A
        temp = deltaUf[:len(fdof)]
        U[fdof, 0] = U[fdof, 0] + temp.T.tolist()[0]
        if (sdof.size > 0):
            U[mdof, 0] = U[mdof, 0] + deltaUf[:(len(fdof))]
            [U[sdof, 0], P_mpc, d2PdUm2_Rs] = MPC_user(U[mdof, 0], F[sdof])

        AllResNrm.append(ResNrm)  # ok<SAGROW>
        AlldUNrm.append(dUNrm)  # ok<SAGROW>

        if (itera == 1) and (iter_load == 0):
            check = checker.build('bigloop.mat')
            check('F', F.T)
            check('ResNrm', ResNrm)
            #check('k_global', k_global)  # KeyError: 'k_global'
            check('Kaa', Kaa)
            check('Pa', Pa.T)
            check('Kff', Kff)
            check('Pf', Pf.T)
            check('deltaUf', deltaUf.T)
            check('U', U.T)

    # Get support reactions

    Fp = F[pdof]

    print ('Load increment ', (iter_load + 1), ' converged after ', itera,
    ' iterations.')
    All_iter.append(itera)  # ok<SAGROW>
    for i in range(nnodes - 1):
        All_soln.append([[U[2 * i, 0]], [U[2 * i + 1, 0]]])



#Compute nodal loads, Von Mises and Tresca
#[StressOut, StressNode] = nodal_stresses(elnodes, stress)
#VonMises = calc_von_mises(StressNode, pois, plane)
#Tresca = calc_tresca(StressNode, pois, plane)

#Write output to text based output file
#write_output_file(file_out, U, displ, Fp, nodes, elnodes, strain, StressOut)

#If GraphOpt=1, start Graphical Output
#if GraphOpt:
#graphical_user_interface(nnodes, coor, nelem, elnodes,
#     StressNode, U, VonMises, Tresca, nloadinc, All_soln)