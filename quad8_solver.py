# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 09:39:32 2014

@author: Dominic
"""
import math
import numpy as np
#import scipy
from solver_utils import read_input_file
from block_diag import block_diag

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

filename = raw_input('Input filename (without extension) ? ')
print filename
[nnodes, ndcoor, nodes, coor, nelem, plane, elnodes, elas, pois, t, ndispl, displ, ncload, cload, nloadinc, mdof, sdof] = read_input_file(filename)
file_out = filename+'.out'

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
    dof[nodes.index(node), dimension-1] = 0
    Up.append(displacement)
Up = np.matrix(Up).T
pdof = np.flatnonzero(dof == 0)
fdof = np.flatnonzero(dof != 0)
#TODO: check if fdof and pdof should actually be a set instead of an array
fdof = np.setdiff1d(fdof[0], mdof)
fdof = np.setdiff1d(fdof[0], sdof)

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
matC = np.asmatrix(block_diag([[c, c * nu], [c * nu, c]], np.eye(2) * c * (1 - nu)))
#print matC

ndnum = range(2,(2 + NodesPerEl))
[colpos, rowpos] = np.meshgrid(range(DofPerEl),range(DofPerEl))

colpos = colpos.flatten()
rowpos = rowpos.flatten()

tol = 3e-5
dUNrm = 1.0

F = np.zeros([2 * nnodes, 1])
LoadFac = (np.array(range(1,nloadinc+1))) / float(nloadinc)
print LoadFac

# FIXME: this 12 should be a variable
stress = np.zeros([nelem, 12])
strain = np.zeros([nelem, 12])


# monster_loop

            # Compute nodal loads, Von Mises and Tresca
#            [StressOut, StressNode] = nodal_stresses(elnodes, stress)
#            VonMises = calc_von_mises(StressNode, pois, plane)
#            Tresca = calc_tresca(StressNode, pois, plane)

            # Write output to text based output file
#            write_output_file(file_out, U, displ, Fp, nodes, elnodes, strain, StressOut)

            # If GraphOpt=1, start Graphical Output
#            if GraphOpt:
#                graphical_user_interface(nnodes, coor, nelem, elnodes, StressNode, U, VonMises, Tresca, nloadinc, All_soln)
#                end