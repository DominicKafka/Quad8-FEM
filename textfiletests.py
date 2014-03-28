from mesh_utils import mesh_output_writer
from mesh_utils import node_coord_displ
from mesh_utils import nodal_forces

#h, l, m, n, option, E, nu, t, load_opt, V0, R, Mag, RadD, filen = mesh_inputs()

h = 1.
l = 10.
m = 10
n = 100
option = 1.
E = 200000.
nu = 0.3
t = 1.
load_opt = 0
Mag = 0
R = 0
RadD = 0
V0 = 100.
filen = 'bob'
DispMat = []
LoadMat = []


coord, displ, elnode, node, el = node_coord_displ(h, l, m, n, load_opt, R)

#forces, forcesb, DispMat, LoadMat = nodal_forces(h, m, V0, load_opt, coord, R, Mag, RadD)

forces, forcesb = nodal_forces(h, m, V0, load_opt, coord, R, Mag, RadD)

mesh_output_writer(filen, node, coord, el, option, elnode, E, nu, t, load_opt, m, DispMat, displ, LoadMat, n, forces, forcesb)


