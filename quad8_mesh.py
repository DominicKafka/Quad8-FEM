# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 12:32:57 2014

@author: Dominic
"""

# clear
# clear(mstring('all'))



from mesh_utils import mesh_inputs
from mesh_utils import node_coord_displ

h, l, m, n, option, E, nu, t, load_opt, V0, R, Mag, RadD, filen = mesh_inputs()

coord, displ, elnode, node, el = node_coord_displ(h, l, m, n, load_opt, R)


#forces, forcesb, DispMat, LoadMat = nodal_forces(h, m, V0, load_opt, coord, R, Mag, RadD)

# DispMat and LoadMat will not be needed for my research project, hence I will 
# confirm whether they need to be in the code at all/whether they need to be translated

forces, forcesb = nodal_forces(h, m, V0, load_opt, coord, R, Mag, RadD)


fid = fopen(filen, mstring('w')); print fid
fprintf(fid, mstring('Number_of_nodes \\n'))
fprintf(fid, mstring('%5d \\n'), node)
fprintf(fid, mstring('Nodal_coordinates \\n'))
fprintf(fid, mstring('%5d %20.15e %20.15e \\n'), coord.cT)
fprintf(fid, mstring('Number_of_elements \\n'))
fprintf(fid, mstring('%5d \\n'), el)
fprintf(fid, mstring('Plane_stress_or_strain \\n'))
fprintf(fid, mstring('%5d \\n'), option)
fprintf(fid, mstring('Element_connectivity \\n'))
fprintf(fid, mstring('%5d %5d %5d %5d %5d %5d %5d %5d %5d \\n'), elnode.cT)
fprintf(fid, mstring('Material_properties \\n'))
fprintf(fid, mstring('%20.15e %20.15e %20.15e \\n'), mcat([E, nu, t]))
fprintf(fid, mstring('Number_of_prescribed_displacements \\n'))
if load_opt == 2:
    fprintf(fid, mstring('%5d \\n'), 2 * node)
elif load_opt != 4:
    fprintf(fid, mstring('%5d \\n'), 4 * m + 2)
else:
    fprintf(fid, mstring('%5d \\n'), 4 * m + 2 + size(DispMat, 1))
end
fprintf(fid, mstring('Prescribed_displacements \\n'))
if load_opt == 2:
    for i in mslice[1:node]:
        fprintf(fid, mstring('%5d 1 %20.15f \\n'), mcat([i, displ(i, 2)]))
        fprintf(fid, mstring('%5d 2 %20.15f \\n'), mcat([i, displ(i, 3)]))
    end
elif load_opt == 3:
    fprintf(fid, mstring('%5d 1 0.0 \\n'), mcat([mslice[1:2 * m + 1]]))
    fprintf(fid, mstring('%5d 2 0.0 \\n'), mcat([mslice[node - (2 * m):node]]))
elif logical_or((load_opt == 1), (load_opt == 0)):
    fprintf(fid, mstring('%5d 1 0.0 \\n'), mcat([mslice[1:2 * m + 1]]))
    fprintf(fid, mstring('%5d 2 0.0 \\n'), mcat([mslice[1:2 * m + 1]]))
elif load_opt == 4:
    fprintf(fid, mstring('%5d 1 0.0 \\n'), mcat([mslice[1:2 * m + 1]]))
    fprintf(fid, mstring('%5d 2 0.0 \\n'), mcat([mslice[node - (2 * m):node]]))
    fprintf(fid, mstring('%5d %5d %20.15f \\n'), DispMat.cT)
end
fprintf(fid, mstring('Number_of_nodal_loads \\n'))
if load_opt == 2:
    fprintf(fid, mstring('0 \\n'))
elif load_opt == 3:
    fprintf(fid, mstring('%5d \\n'), size(LoadMat, 1))
elif logical_or((load_opt == 1), (load_opt == 0)):
    fprintf(fid, mstring('%5d \\n'), 2 * m + 1)
elif load_opt == 4:
    fprintf(fid, mstring(' 0 \\n'))
end
fprintf(fid, mstring('Nodal_loads \\n'))
if load_opt == 0:
    for i in mslice[1:2 * m + 1]:
        c_node = (3 * m + 2) * n + i
        fprintf(fid, mstring('%5d 2 %20.15f \\n'), mcat([c_node, forces(i)]))
    end
elif load_opt == 1:
    for i in mslice[1:2 * m + 1]:
        c_node = (3 * m + 2) * n + i
        fprintf(fid, mstring('%5d 1 %20.15f \\n'), mcat([c_node, forcesb(i)]))
    end
elif load_opt == 3:
    fprintf(fid, mstring('%5d %5d %20.15f \\n'), LoadMat.cT)
end
fprintf(fid, mstring('Number_of_load_increments \\n'))
fprintf(fid, mstring(' 10 \\n'))
fprintf(fid, mstring('Number_of_MPCs \\n'))
fprintf(fid, mstring('  0 \\n'))
fclose(fid)

figure(1)
plot(coord(mslice[:], 2), coord(mslice[:], 3), mstring('x'))
for i in mslice[1:node]:
    text(coord(i, 2), coord(i, 3), num2str(i))
end
axis(mstring('equal'))
