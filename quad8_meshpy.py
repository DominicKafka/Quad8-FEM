# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 12:32:57 2014

@author: Dominic
"""

# clear
# clear(mstring('all'))


import math

h = float(input('Height of rectangle ? '))
l = float(input('Length of rectangle ? '))
m = int(input('Increments along height ? '))
n = int(input('Increments along length ? '))
option = int(input('0:Plane strain or 1:Plane stress ? '))
E = int(input('Elasticity modulus ? '))
nu = float(input('Poisson\'s ratio ? '))
t = float(input('Element thickness ? '))
print('0: Parabolic shear stress or')
print('1: Linear bending stress at beam tip or')
print('2: Beam forced into radius or')
print('3: Beam transformed to 90 degree section loaded with pressure')
load_opt = int(input('4: Beam transformed to 90 degree section, inside surface displaced '))

if load_opt == 0:
    V0 = float(input('Magnitude of total shear force at beam tip ? '))
elif load_opt == 1:
    V0 = float(input('Magnitude of bending stress at top surface at beam tip ? '))
elif load_opt == 2:
    R = float(input('Enter radius that beam is deformed into. '))
elif load_opt == 3:
    R = 2. * l / math.pi
    Mag = float(input('Internal pressure (in undeformed configuration) ? '))
elif load_opt == 4:
    R = 2 * l / math.pi
    RadD = float(input('Radial displacement ? '))

filen = raw_input('Write output to which file ? ')

coord = []
displ = []
node = 0
deltx = l / (2 * n)
for i in range(2 * n + 1):
    if (((i+1) % 2) == 1):
        delty = h / (2 * m)
        for j in range(2 * m + 1):
            node = node + 1
            x = i * deltx
            y = j * delty - h / 2
            coord.append([node,x,y])
            if ((load_opt == 2) or (load_opt == 3) or (load_opt == 4)):
                u = (R - y) * math.sin(x / R) - x
                v = R - (R - y) * math.cos(x / R) - y
                if load_opt == 2:
                    displ.append([node, u, v])
                    coord.append([node, x, y])
                else:
                    coord.append([node(x + u), (y + v)])

    else:
        delty = h / m
        for j in range(m + 1):
            node = node + 1
            x = (i - 1) * deltx
            y = (j - 1) * delty - h / 2
            coord.append([node,x,y])
            if ((load_opt == 2) or (load_opt == 3) or (load_opt == 4)):
                u = (R - y) * math.sin(x / R) - x
                v = R - (R - y) * math.cos(x / R) - y
                if load_opt == 2:
                    displ.append([node, u, v])
                    coord.append([node, x, y])
                else:
                    coord.append([node(x + u), (y + v)])

el = 0
for i in mslice[1:n]:
    for j in mslice[1:m]:
        elnode1 = (3 * m + 2) * (i - 1) + 2 * (j - 1) + 1
        elnode2 = (3 * m + 2) * (i - 1) + 2 * (j - 1) + 3
        elnode3 = (3 * m + 2) * i + 2 * (j - 1) + 1
        elnode4 = (3 * m + 2) * i + 2 * (j - 1) + 3
        elnode5 = 2 * m + 1 + (3 * m + 2) * (i - 1) + j
        elnode6 = (3 * m + 2) * i + 2 * (j - 1) + 2
        elnode7 = 2 * m + 1 + (3 * m + 2) * (i - 1) + j + 1
        elnode8 = (3 * m + 2) * (i - 1) + 2 * (j - 1) + 2
        el1 = el + 1
        elnode(el1, mslice[:]).lvalue = mcat([el1, elnode1, elnode3, elnode4, elnode2, elnode5, elnode6, elnode7, elnode8])
        el = el + 1
    end
end

#Compute equivalant nodal forces applied at beam tip
delty = h / m
forces = zeros(2 * m + 1, 1)
forcesb = zeros(2 * m + 1, 1)
if (load_opt < 2):
    for i in mslice[1:m]:
        y1 = -h / 2 + (i - 1) * delty
        y3 = -h / 2 + i * delty
        y2 = -h / 2 + (i - 0.5) * delty
        # Compute equivalent nodal loads for parabolic shear stress
        FA = -1 / 3 * (5 * y1 ** 3 - 11 * y2 * y1 ** 2 + 3 * y3 * y1 ** 2 + 6 * y1 * y2 ** 2 - 8 * y1 * y2 * y3 + 3 * y1 * y3 ** 2 + y3 ** 3 - 5 * y2 * y3 ** 2 + 6 * y3 * y2 ** 2) * V0 / (-y3 + y1) / h
        FB = 1 / 3 * (3 * y1 ** 3 + 9 * y3 * y1 ** 2 - 16 * y2 * y1 ** 2 + 12 * y1 * y2 ** 2 + 9 * y1 * y3 ** 2 - 16 * y1 * y2 * y3 + 3 * y3 ** 3 + 12 * y3 * y2 ** 2 - 16 * y2 * y3 ** 2) * V0 / (-y3 + y1) / h
        FC = -1 / 3 * (y1 ** 3 - 5 * y2 * y1 ** 2 + 3 * y3 * y1 ** 2 + 6 * y1 * y2 ** 2 - 8 * y1 * y2 * y3 + 3 * y1 * y3 ** 2 + 5 * y3 ** 3 - 11 * y2 * y3 ** 2 + 6 * y3 * y2 ** 2) * V0 / (-y3 + y1) / h
        # Compute equivalant nodal loads for linear bending stress
        F1 = 1 / 20 * (-18 * y3 ** 4 - 48 * y1 * y3 ** 3 + 80 * y3 ** 3 * y2 - 80 * y3 ** 2 * y2 ** 2 + 5 * h ** 2 * y3 ** 2 + 120 * y3 ** 2 * y2 * y1 - 48 * y3 ** 2 * y1 ** 2 - 30 * y3 * y2 * h ** 2 + 20 * h ** 2 * y3 * y1 - 80 * y3 * y1 * y2 ** 2 + 120 * y3 * y2 * y1 ** 2 - 48 * y3 * y1 ** 3 + 60 * y2 ** 2 * h ** 2 - 78 * y1 ** 4 + 35 * h ** 2 * y1 ** 2 - 90 * y2 * h ** 2 * y1 + 160 * y2 * y1 ** 3 - 80 * y1 ** 2 * y2 ** 2) * V0 / (y3 - y1) / h ** 3
        F2 = -1 / 10 * (-28 * y3 ** 4 - 68 * y1 * y3 ** 3 + 120 * y3 ** 3 * y2 - 80 * y3 ** 2 * y2 ** 2 + 5 * h ** 2 * y3 ** 2 + 120 * y3 ** 2 * y2 * y1 - 48 * y3 ** 2 * y1 ** 2 - 68 * y3 * y1 ** 3 + 50 * h ** 2 * y3 * y1 - 60 * y3 * y2 * h ** 2 - 80 * y3 * y1 * y2 ** 2 + 120 * y3 * y2 * y1 ** 2 + 5 * h ** 2 * y1 ** 2 - 28 * y1 ** 4 + 60 * y2 ** 2 * h ** 2 - 60 * y2 * h ** 2 * y1 - 80 * y1 ** 2 * y2 ** 2 + 120 * y2 * y1 ** 3) * V0 / (y3 - y1) / h ** 3
        F3 = 1 / 20 * (-78 * y3 ** 4 - 48 * y1 * y3 ** 3 + 160 * y3 ** 3 * y2 - 80 * y3 ** 2 * y2 ** 2 + 35 * h ** 2 * y3 ** 2 + 120 * y3 ** 2 * y2 * y1 - 48 * y3 ** 2 * y1 ** 2 - 90 * y3 * y2 * h ** 2 + 20 * h ** 2 * y3 * y1 - 80 * y3 * y1 * y2 ** 2 + 120 * y3 * y2 * y1 ** 2 - 48 * y3 * y1 ** 3 + 60 * y2 ** 2 * h ** 2 - 18 * y1 ** 4 + 5 * h ** 2 * y1 ** 2 - 30 * y2 * h ** 2 * y1 + 80 * y2 * y1 ** 3 - 80 * y1 ** 2 * y2 ** 2) * V0 / (y3 - y1) / h ** 3
        # forces contains the equivalent nodal loads for
        # parabolic shear stress distribution at beam tip
        forces(mslice[2 * (i - 1) + 1:2 * i + 1]).lvalue = forces(mslice[2 * (i - 1) + 1:2 * i + 1]) + mcat([F1, F2, F3]).cT
        # forcesb contains the equivalant nodal loads for
        # linear bending stress applied at the beam tip
        forcesb(mslice[2 * (i - 1) + 1:2 * i + 1]).lvalue = forcesb(mslice[2 * (i - 1) + 1:2 * i + 1]) + mcat([FA, FB, FC]).cT
    end
end

if load_opt == 3:
    Interior = find(sqrt(coord(mslice[:], 2) **elpow** 2 + (coord(mslice[:], 3) - R) **elpow** 2) - (R - h / 2) < 0.001)
    LoadMat = mcat([Interior(1), 2 - Mag / 6])
    for i in mslice[2:(length(Interior) - 1)]:
        dx = coord(Interior(i), 2)
        dy = coord(Interior(i), 3) - R
        dL = sqrt(dx ** 2 + dy ** 2)
        if mod(i, 2) == 0:
            Loadx = 2 / 3 * Mag * dx / dL
            Loady = 2 / 3 * Mag * dy / dL
        else:
            Loadx = 1 / 3 * Mag * dx / dL
            Loady = 1 / 3 * Mag * dy / dL
        end

        Interior(i)
        1

        Interior(i)
        2
    end

    Interior(end)
    1
end

if load_opt == 4:
    Interior = find(sqrt(coord(mslice[:], 2) **elpow** 2 + (coord(mslice[:], 3) - R) **elpow** 2) - (R - h / 2) < 0.001)
    Exterior = find(sqrt(coord(mslice[:], 2) **elpow** 2 + (coord(mslice[:], 3) - R) **elpow** 2) > (R + h / 2) - 0.001)

    DispMat = mcat([Interior(1), 2 - RadD])
    for i in mslice[2:(length(Interior) - 1)]:
        dx = coord(Interior(i), 2)
        dy = coord(Interior(i), 3) - R
        dL = sqrt(dx ** 2 + dy ** 2)
        Dx = RadD * dx / dL
        Dy = RadD * dy / dL

        Interior(i)
        1
        Dx
        Interior(i)
        2
    end

    Interior(end)
    1
end

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