function gen_mesh_output(filen,node,coord,option,elnode,E,nu,t,load_opt,m,n,DispMat,l,displ,LoadMat,forces,forcesb)

fid=fopen(filen,'w');
fprintf(fid,'Number_of_nodes \n');
fprintf(fid,'%5d \n',node);
fprintf(fid,'Nodal_coordinates \n');
fprintf(fid,'%5d %20.15e %20.15e \n',coord');
fprintf(fid,'Number_of_elements \n');
fprintf(fid,'%5d \n',el);
fprintf(fid,'Plane_stress_or_strain \n');
fprintf(fid,'%5d \n',option);
fprintf(fid,'Element_connectivity \n');
fprintf(fid,'%5d %5d %5d %5d %5d %5d %5d %5d %5d \n',elnode');
fprintf(fid,'Material_properties \n');
fprintf(fid,'%20.15e %20.15e %20.15e \n',[E nu t]);
fprintf(fid,'Number_of_prescribed_displacements \n');
if load_opt==2
    fprintf(fid,'%5d \n',2*node);
elseif load_opt~=4
    fprintf(fid,'%5d \n',4*m+2);
else
    fprintf(fid,'%5d \n',4*m+2+size(DispMat,1));
end
fprintf(fid,'Prescribed_displacements \n');
if load_opt==2
    for i=1:node
        fprintf(fid,'%5d 1 %20.15f \n',[i displ(i,2)]);
        fprintf(fid,'%5d 2 %20.15f \n',[i displ(i,3)]);
    end
elseif load_opt==3
    fprintf(fid,'%5d 1 0.0 \n',1:2*m+1);
    fprintf(fid,'%5d 2 0.0 \n',node-(2*m):node);
elseif any(load_opt == [0, 1])
    fprintf(fid,'%5d 1 0.0 \n',1:2*m+1);
    fprintf(fid,'%5d 2 0.0 \n',1:2*m+1);
elseif load_opt==4
    fprintf(fid,'%5d 1 0.0 \n',1:2*m+1);
    fprintf(fid,'%5d 2 0.0 \n',node-(2*m):node);
    fprintf(fid,'%5d %5d %20.15f \n',DispMat');
end
fprintf(fid,'Number_of_nodal_loads \n');
if load_opt==2
    fprintf(fid,'0 \n');
elseif load_opt==3
    fprintf(fid,'%5d \n',size(LoadMat,1));
elseif any(load_opt == [0, 1])
    fprintf(fid,'%5d \n',2*m+1);
elseif load_opt==4
    fprintf(fid,' 0 \n');
end;
fprintf(fid,'Nodal_loads \n');
if load_opt==0
    for i=1:2*m+1
        c_node = (3*m+2)*n+i;
        fprintf(fid,'%5d 2 %20.15f \n',[c_node forces(i)]);
    end
elseif load_opt==1
    for i=1:2*m+1
        c_node = (3*m+2)*n+i;
        fprintf(fid,'%5d 1 %20.15f \n',[c_node forcesb(i)]);
    end
elseif load_opt==3
    fprintf(fid,'%5d %5d %20.15f \n',LoadMat');
end
fprintf(fid,'Number_of_load_increments \n');
fprintf(fid,' 10 \n');
fprintf(fid,'Number_of_MPCs \n');
fprintf(fid,'  0 \n');
fclose(fid);

figure(1)
plot(coord(:,2),coord(:,3),'x')
for i=1:node
    text(coord(i,2),coord(i,3),num2str(i));
end
axis equal
