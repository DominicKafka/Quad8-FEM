# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 14:01:10 2014

@author: Dominic
"""

def readsectionfile(filename):
    """ Read a file containing string section names and array-like rows of numbers

    This returns a dictionary with an entry for each section and parses the array entries into floats as a list of rows.
    """
    section = {}
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

    import time
    from numpy import array

    tic = time.time()

    section = readsectionfile(filename)

    # Read number of nodes
    [[nnodes]] = section['Number_of_nodes']

    #Read nodal coordinates
    ndcoor = array(section['Nodal_coordinates'])
    nodes = ndcoor[:, 0]
    coor = ndcoor[:, 1:]

    print nodes
    print coor

    # Read number of elements
    [[nelem]] = section['Number_of_elements']

    # Read if elements are plane stress or plain strain
    [[plane]] = section['Plane_stress_or_strain']

    # Read element number and element connectivity
    elnodes = array(section['Element_connectivity'], int)

    print nelem
    print plane
    print elnodes

    # Read material constants
    [[elas, pois, t]] = section['Material_properties']
    
    print elas
    print pois
    print t    

    # Read number of prescribed displacements
    [[ndispl]] = section['Number_of_prescribed_displacements']

    # Read prescribed displacements
    displ = array(section['Prescribed_displacements'])
    print displ
    
    # Read number of  nodal loads
    [[ncload]] = section['Number_of_nodal_loads']
    print ncload

    # Read nodal loads
    cload = array(section['Nodal_loads'])
    print cload
    
    # Read number of load increments
    [[nloadinc]] = section['Number_of_load_increments']
    print nloadinc

   
    # Read number of MPCs
    [[nMPC]] = section['Number_of_MPCs']
    print nMPC
    

    MasterDOF = []
    SlaveDOF = []
    
    toc = time.time()
    time = toc-tic
    print time    
    print 'Done reading input file '+filename+ ' in ' +str(time)+' seconds'
    return nnodes, ndcoor, nodes, coor, nelem, plane, elnodes, elas, pois, t, ndispl, displ, ncload, cload, nloadinc, MasterDOF, SlaveDOF
