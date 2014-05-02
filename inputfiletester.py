# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 18:45:58 2014

@author: Dominic
"""

from solver_utils import read_input_file

filename = "Beam2by20.inp"

(ndcoor, nodes, coor, plane, elnodes, elas, pois, t, 
     displ, cload, nloadinc, MasterDOF, SlaveDOF) = (read_input_file(filename))
