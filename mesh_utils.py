def mesh_inputs():
    
    import math    
    V0 = 0
    R = 0
    Mag = 0
    RadD = 0

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
    return h, l, m, n, option, E, nu, t, load_opt, V0, R, Mag, RadD, filen
    
    
def node_coord_displ(h, l, m, n, load_opt, R):

    import math
    coord = []
    displ = []
    elnode = []
    node = 0
    deltx = l / (2. * n)
    for i in range(2 * n + 1):
        multiplier = ((i + 1) % 2) + 1
        delty = h / (float(multiplier) * m)
        for j in range(multiplier * m + 1):
            node = node + 1
            x = i * deltx
            y = j * delty - h / 2.
            coord.append([node, x, y])
            if load_opt in [2, 3, 4]:
                u = (R - y) * math.sin(x / R) - x
                v = R - (R - y) * math.cos(x / R) - y
                if load_opt == 2:
                    displ.append([node, u, v])
                    coord.append([node, x, y])
                else:
                    coord.append([node(x + u), (y + v)])

    el = 0
    for i in range(n):
        for j in range(m):
            elnode1 = (3 * m + 2) * (i) + 2 * (j) + 1
            elnode2 = (3 * m + 2) * (i) + 2 * (j) + 3
            elnode3 = (3 * m + 2) * (i + 1) + 2 * (j) + 1
            elnode4 = (3 * m + 2) * (i + 1) + 2 * (j) + 3
            elnode5 = 2 * m + 1 + (3 * m + 2) * (i) + (j + 1)
            elnode6 = (3 * m + 2) * (i + 1) + 2 * (j) + 2
            elnode7 = 2 * m + 1 + (3 * m + 2) * (i) + (j + 1) + 1
            elnode8 = (3 * m + 2) * (i) + 2 * (j) + 2
            el1 = el + 1
            elnode.append([el1, elnode1, elnode3, elnode4, elnode2, elnode5, elnode6, elnode7, elnode8])
            el = el + 1

    return coord, displ, elnode, node
    

def nodal_forces(h, m, V0, load_opt, coord, R, Mag, RadD):
        
    import numpy
    import math
    DispMat = []
    LoadMat = []    
    
    #Compute equivalant nodal forces applied at beam tip
    delty = h / float(m)
    forces = numpy.array(numpy.zeros(shape=(2 * m + 1, 1)), dtype=float)
    forcesb = numpy.array(numpy.zeros(shape=(2 * m + 1, 1)), dtype=float)
    if (load_opt < 2):
        for i in range(m):
            y1 = -h / 2. + (i) * delty
            y3 = -h / 2. + (i + 1) * delty
            y2 = -h / 2. + (i + 0.5) * delty
            # Compute equivalent nodal loads for parabolic shear stress
            FA = -1. / 3 * (5 * y1 ** 3 - 11 * y2 * y1 ** 2 + 3 * y3 * y1 ** 2 + 6 * y1 * y2 ** 2 - 8 * y1 * y2 * y3 + 3 * y1 * y3 ** 2 + y3 ** 3 - 5 * y2 * y3 ** 2 + 6 * y3 * y2 ** 2) * V0 / (-y3 + y1) / h
            FB = 1. / 3 * (3 * y1 ** 3 + 9 * y3 * y1 ** 2 - 16 * y2 * y1 ** 2 + 12 * y1 * y2 ** 2 + 9 * y1 * y3 ** 2 - 16 * y1 * y2 * y3 + 3 * y3 ** 3 + 12 * y3 * y2 ** 2 - 16 * y2 * y3 ** 2) * V0 / (-y3 + y1) / h
            FC = -1. / 3 * (y1 ** 3 - 5 * y2 * y1 ** 2 + 3 * y3 * y1 ** 2 + 6 * y1 * y2 ** 2 - 8 * y1 * y2 * y3 + 3 * y1 * y3 ** 2 + 5 * y3 ** 3 - 11 * y2 * y3 ** 2 + 6 * y3 * y2 ** 2) * V0 / (-y3 + y1) / h
            # Compute equivalant nodal loads for linear bending stress
            F1 = 1. / 20 * (-18 * y3 ** 4 - 48 * y1 * y3 ** 3 + 80 * y3 ** 3 * y2 - 80 * y3 ** 2 * y2 ** 2 + 5 * h ** 2 * y3 ** 2 + 120 * y3 ** 2 * y2 * y1 - 48 * y3 ** 2 * y1 ** 2 - 30 * y3 * y2 * h ** 2 + 20 * h ** 2 * y3 * y1 - 80 * y3 * y1 * y2 ** 2 + 120 * y3 * y2 * y1 ** 2 - 48 * y3 * y1 ** 3 + 60 * y2 ** 2 * h ** 2 - 78 * y1 ** 4 + 35 * h ** 2 * y1 ** 2 - 90 * y2 * h ** 2 * y1 + 160 * y2 * y1 ** 3 - 80 * y1 ** 2 * y2 ** 2) * V0 / (y3 - y1) / h ** 3
            F2 = -1. / 10 * (-28 * y3 ** 4 - 68 * y1 * y3 ** 3 + 120 * y3 ** 3 * y2 - 80 * y3 ** 2 * y2 ** 2 + 5 * h ** 2 * y3 ** 2 + 120 * y3 ** 2 * y2 * y1 - 48 * y3 ** 2 * y1 ** 2 - 68 * y3 * y1 ** 3 + 50 * h ** 2 * y3 * y1 - 60 * y3 * y2 * h ** 2 - 80 * y3 * y1 * y2 ** 2 + 120 * y3 * y2 * y1 ** 2 + 5 * h ** 2 * y1 ** 2 - 28 * y1 ** 4 + 60 * y2 ** 2 * h ** 2 - 60 * y2 * h ** 2 * y1 - 80 * y1 ** 2 * y2 ** 2 + 120 * y2 * y1 ** 3) * V0 / (y3 - y1) / h ** 3
            F3 = 1. / 20 * (-78 * y3 ** 4 - 48 * y1 * y3 ** 3 + 160 * y3 ** 3 * y2 - 80 * y3 ** 2 * y2 ** 2 + 35 * h ** 2 * y3 ** 2 + 120 * y3 ** 2 * y2 * y1 - 48 * y3 ** 2 * y1 ** 2 - 90 * y3 * y2 * h ** 2 + 20 * h ** 2 * y3 * y1 - 80 * y3 * y1 * y2 ** 2 + 120 * y3 * y2 * y1 ** 2 - 48 * y3 * y1 ** 3 + 60 * y2 ** 2 * h ** 2 - 18 * y1 ** 4 + 5 * h ** 2 * y1 ** 2 - 30 * y2 * h ** 2 * y1 + 80 * y2 * y1 ** 3 - 80 * y1 ** 2 * y2 ** 2) * V0 / (y3 - y1) / h ** 3
            # forces contains the equivalent nodal loads for
            # parabolic shear stress distribution at beam tip
            forces[2 * (i),0] = (forces[2 * (i),0]) + F1
            forces[2 * (i) + 1,0] = (forces[2 * (i) + 1,0]) + F2
            forces[2 * (i) + 2,0] = (forces[2 * (i) + 2,0]) + F3
            # forcesb contains the equivalant nodal loads for
            # linear bending stress applied at the beam tip
            forcesb[2 * (i),0] = (forcesb[2 * (i),0]) + FA
            forcesb[2 * (i) + 1,0] = (forcesb[2 * (i) + 1,0]) + FB
            forcesb[2 * (i) + 2,0] = (forcesb[2 * (i) + 2,0]) + FC
            
            
    if load_opt == 3:
        Interior = numpy.where(math.sqrt(coord[:, 2]** 2 + (coord[:, 3] - R)** 2) - (R - h / 2.) < 0.001)
        LoadMat = [Interior[1], 2, -Mag / 6]
        for i in range(1,(len(Interior))):
            dx = coord[Interior[i], 1]
            dy = coord[Interior[i], 2] - R
            dL = math.sqrt(dx ** 2 + dy ** 2)
            if (i % 2) == 0:
                Loadx = 2. / 3 * Mag * dx / dL
                Loady = 2. / 3 * Mag * dy / dL
            else:
                Loadx = 1. / 3 * Mag * dx / dL
                Loady = 1. / 3 * Mag * dy / dL
                
            LoadMat = [LoadMat, Interior(i), 1, Loadx]
            LoadMat = [LoadMat, Interior(i), 2, Loady]
                
        LoadMat = [LoadMat, Interior.end, 1, 1./6*Mag];
    
    if load_opt == 4:
        Interior = numpy.where(math.sqrt(coord[:, 2]** 2 + (coord[:, 3] - R)** 2) - (R - h / 2.) < 0.001)
        Exterior = numpy.where(math.sqrt(coord[:, 2]** 2 + (coord[:, 3] - R)** 2) > (R + h / 2.) - 0.001)
        
        DispMat = ([Interior(1), 2, -RadD])
        for i in range(1,(len(Interior))):
            dx = coord[Interior[i], 1]
            dy = coord[Interior[i], 2] - R
            dL = math.sqrt(dx ** 2 + dy ** 2)
            Dx = RadD * dx / dL
            Dy = RadD * dy / dL
            
            DispMat = [DispMat, Interior[i], 1, Dx, Interior[i], 2, Dy]
        DispMat = [DispMat, Interior.end, 1, RadD]    
        
    return forces, forcesb, DispMat, LoadMat
#    return forces, forcesb


def mesh_output_writer(filen, node, coord, el, option, elnode, E, nu, t, load_opt, m, DispMat, displ, LoadMat, n, forces, forcesb):
    
    import pylab as pl
    
    fid = open(filen,'w')
    fid.write("Number_of_nodes \n")
    fid.write("%5f \n"% node)
    fid.write("Nodal_coordinates \n")
    for i in range(node):
        fid.write("%f %f %f \n"% tuple(coord[i]))
    fid.write('Number_of_elements \n')
    fid.write('%5d \n'% el)
    fid.write('Plane_stress_or_strain \n')
    fid.write('%5d \n' % option)
    fid.write('Element_connectivity \n')
    for j in range(el):    
        fid.write('%5d %5d %5d %5d %5d %5d %5d %5d %5d \n'% tuple(elnode[j]))
    fid.write('Material_properties \n')
    fid.write('%f %f %f \n'% tuple([E, nu, t]))
    fid.write('Number_of_prescribed_displacements \n')
    if load_opt == 2:
        fid.write('%5d \n'% 2 * node)
    elif load_opt == 4:
        fid.write('%5d \n'% 4 * m + 2)
    else:
        fid.write('%5d \n'% (4 * m + 2 + int(len(DispMat))))
    fid.write('Prescribed_displacements \n')
    if load_opt == 2:
        for i in range(node):
            fid.write('%5d 1 %20.15f \n'% (i, displ[i, 1]))
            fid.write('%5d 2 %20.15f \n'% (i, displ[i, 2]))

    elif load_opt ==  3:
        for i in range(2 * m + 1):
            fid.write('%5d 1 0.0 \n'% (i+1))
        for j in range((node - (2 * m)),node):
            fid.write('%5d 2 0.0 \n'% (j+1))
    elif load_opt in [0, 1]:
        for i in range(2 * m + 1):
            fid.write('%5d 1 0.0 \n'% int(i+1))
        for j in range(2 * m + 1):
            fid.write('%5d 2 0.0 \n'% int(j+1))
    elif load_opt == 4:
        for i in range(2 * m + 1):        
            fid.write('%5d 1 0.0 \n'% int(i+1))
        for j in range((node - (2 * m)),node):
            fid.write('%5d 2 0.0 \n'% int(j+1))
        for k in range(int(len(DispMat))):
            fid.write('%5d %5d %f \n'% DispMat[k])

    fid.write('Number_of_nodal_loads \n')
    if load_opt == 2:
        fid.write('0 \n')
    elif load_opt == 3:
        fid.write('%5d \n' % int(len(LoadMat)))
    elif load_opt in [0, 1]:
        fid.write('%5d \n'% (2 * m + 1))
    elif load_opt == 4:
        fid.write(' 0 \n')
    fid.write('Nodal_loads \n')
                                                    
    if load_opt == 0:
        for i in range(2 * m + 1):
            c_node = (3 * m + 2) * n + (i+1)
            fid.write('%5d 2 %20.15f \n'% tuple([c_node, forces[i]]))

    elif load_opt == 1:
        for i in range(2 * m + 1):
            c_node = (3 * m + 2) * n + (i+1)
            fid.write('%5d 1 %20.15f \n'% tuple([c_node, forcesb[i]]))
            
    elif load_opt == 3:
        for i in range(len(LoadMat)):        
            fid.write('%5d %5d %20.15f \n'% LoadMat[i])
   
    fid.write('Number_of_load_increments \n')
    fid.write(' 10 \n')
    fid.write('Number_of_MPCs \n')
    fid.write('  0 \n')
    fid.close()
                                                                        
    pl.figure(1)
    label = [x[0] for x in coord]    
    x_val = [x[1] for x in coord]
    y_val = [x[2] for x in coord]
    pl.plot(x_val, y_val,'x')
    for i in range(node):
        pl.text(x_val[i], y_val[i], label[i])
    pl.axis('equal')
    pl.show()
    
    
    
    
    
    