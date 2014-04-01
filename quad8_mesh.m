clear all

[h,l,m,n,option,E,nu,t,load_opt,V0,R,Mag,RadD,filen] = mesh_inputs;

[coord,displ,elnode,node] = node_coord_displ(h,l,m,n,load_opt);

[forces,forcesb,DispMat,LoadMat] = nodal_forces(h,m,V0,load_opt,coord,R,Mag,RadD); 

gen_mesh_output(filen,node,coord,option,elnode,E,nu,t,load_opt,m,n,DispMat,l,displ,LoadMat,forces,forcesb)

