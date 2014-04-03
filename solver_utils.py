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
    """ Read sectioned input file for a FEM problem"""

    #TODO: the proper types should be used rather than just float or int for everything
    import time
    from numpy import array

    tic = time.time()

    section = readsectionfile(filename)

    scalarsections = ['Number_of_nodes',
                      'Number_of_elements',
                      'Plane_stress_or_strain',
                      'Number_of_prescribed_displacements',
                      'Number_of_nodal_loads',
                      'Number_of_load_increments',
                      'Number_of_MPCs',
                      ]

    nnodes, nelem, plane, ndispl, ncload, nloadinc, nMPC = [int(section[s][0][0]) for s in scalarsections]
    # TODO: plane should be of type bool

    arraysections = ['Nodal_coordinates',
                     'Element_connectivity',
                     'Nodal_loads',
                     ]    

    ndcoor, elnodes, cload = [array(section[s]) for s in arraysections]
    # TODO: elnodes should be integer

    # sort the rows of displ - keep it as a list
    displ = sorted(section['Prescribed_displacements'])

    nodes = ndcoor[:, 0]
    coor = ndcoor[:, 1:]

    [[elas, pois, t]] = section['Material_properties']

    MasterDOF = []
    SlaveDOF = []
    
    toc = time.time()
    time = toc-tic
    print 'Done reading input file {} in {} seconds'.format(filename, time)
    return nnodes, ndcoor, nodes, coor, nelem, plane, elnodes, elas, pois, t, ndispl, displ, ncload, cload, nloadinc, MasterDOF, SlaveDOF
