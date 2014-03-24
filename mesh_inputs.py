def mesh_inputs():
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
    R = 0

    filen = raw_input('Write output to which file ? ')
    return h, l, m, n, option, E, nu, t, load_opt, V0, R, Mag, RadD, filen
