# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 12:32:57 2014

@author: Dominic
"""

from mesh_utils import mesh_inputs
from mesh_utils import mesh_output_writer
from mesh_utils import node_coord_displ
from mesh_utils import nodal_forces

h, l, m, n, option, E, nu, t, load_opt, V0, R, Mag, RadD, filen = mesh_inputs()

coord, displ, elnode, node, el = node_coord_displ(h, l, m, n, load_opt, R)

forces, forcesb, DispMat, LoadMat = nodal_forces(h, m, V0, load_opt, coord, R, Mag, RadD)

mesh_output_writer(filen, node, coord, el, option, elnode, E, nu, t, load_opt, m, DispMat, displ, LoadMat, n, forces, forcesb)