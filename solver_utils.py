# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 14:01:10 2014

@author: Dominic
"""

def read_input_file(filename, NodesPerEl):
    
    import time
    from numpy import array
    # Open input file
    fid = open(filename, 'r')

    tic = time.time()
    # Read number of nodes
    for i in range(2):    
        line = fid.readline()
        nnodes = line.split()    
    nnodes = int(float(nnodes[0]))
    print nnodes
    
    #Read nodal coordinates
    ndcoor = []
    line = fid.readline()  
    for i in range(nnodes):
        line = fid.readline()        
        vector_list = line.split()
        ndcoor.append(vector_list)
    # convert the list to an array
    ndcoor = array(ndcoor, dtype=float) 
    
    print ndcoor
    
    nodes = ndcoor[:,0]
    coor = ndcoor[:,1:3]
    
    print nodes
    print coor
    

    # Read number of elements
    for i in range(2):   
        line = fid.readline()
    nelem = line.split()
    nelem = int(nelem[0])

    # Read if elements are plane stress or plain strain
    for i in range(2):   
        line = fid.readline()
    plane = line.split()
    plane = int(plane[0])
        
    # Read element number and element connectivity
    line = fid.readline()
    elnodes = []
    for i in range(nelem):        
        line = fid.readline()        
        vector_list = line.split()
        elnodes.append(vector_list)
    # convert the list to an array
    elnodes = array(elnodes, dtype=float)
    
    print nelem
    print plane
    print elnodes
    
    # Read material constants
    for i in range(2):   
        line = fid.readline()
    elas, nu, t = line.split()
    elas = float(elas)
    pois = float(nu)
    t = float(t)
    
    print elas
    print pois
    print t    

    # Read number of prescribed displacements
    for i in range(2):   
        line = fid.readline()
    ndispl = line.split()
    ndispl = int(ndispl[0])
    print ndispl

    # Read prescribed displacements
    line = fid.readline()
    displ = []
    for i in range(ndispl):        
        line = fid.readline()        
        vector_list = line.split()
        displ.append(vector_list)
    # convert the list to an array
    displ = array(displ, dtype=float)
    print displ
    
    # Read number of  nodal loads
    for i in range(2):   
        line = fid.readline()
    ncload = line.split()
    ncload = int(ncload[0])
    print ncload

    # Read nodal loads
    line = fid.readline()
    cload = []
    for i in range(ncload):        
        line = fid.readline()        
        vector_list = line.split()
        cload.append(vector_list)
    # convert the list to an array
    cload = array(cload, dtype=float)
    print cload
    
    # Read number of load increments
    for i in range(2):   
        line = fid.readline()
    nloadinc = line.split()
    nloadinc = int(nloadinc[0])
    print nloadinc

   
    # Read number of MPCs
    for i in range(2):   
        line = fid.readline()
    nMPC = line.split()
    nMPC = int(nMPC[0])
    print nMPC
    

    MasterDOF = []
    SlaveDOF = []
    
    fid.close()
    toc = time.time()
    time = toc-tic
    print time    
    print 'Done reading input file '+filename+ ' in ' +str(time)+' seconds'
    return nnodes, ndcoor, nodes, coor, nelem, plane, elnodes, elas, pois, t, ndispl, displ, ncload, cload, nloadinc, MasterDOF, SlaveDOF
