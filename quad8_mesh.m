clear all

[h,l,m,n,option,E,nu,t,load_opt,V0,R,Mag,RadD,filen] = mesh_inputs

[coord,displ,elnode,node,el] = node_coord_displ(h,l,m,n,load_opt);


%Compute equivalant nodal forces applied at beam tip
delty = h / m;
forces = zeros(2*m+1,1);
forcesb = zeros(2*m+1,1);
if (load_opt<2)
    for i=1:m;
        y1 = -h/2 + (i-1)*delty;
        y3 = -h/2 + i*delty;
        y2 = -h/2 + (i-0.5)*delty;
        % Compute equivalent nodal loads for parabolic shear stress
        FA = -1/3*(5*y1^3-11*y2*y1^2+3*y3*y1^2+6*y1*y2^2-8*y1*y2*y3+3*y1*y3^2+ ...
            y3^3-5*y2*y3^2+6*y3*y2^2)*V0/(-y3+y1)/h;
        FB =  1/3*(3*y1^3+9*y3*y1^2-16*y2*y1^2+12*y1*y2^2+9*y1*y3^2-16*y1*y2*y3+ ...
            3*y3^3+12*y3*y2^2-16*y2*y3^2)*V0/(-y3+y1)/h;
        FC = -1/3*(y1^3-5*y2*y1^2+3*y3*y1^2+6*y1*y2^2-8*y1*y2*y3+3*y1*y3^2+ ...
            5*y3^3-11*y2*y3^2+6*y3*y2^2)*V0/(-y3+y1)/h;
        % Compute equivalant nodal loads for linear bending stress
        F1 = 1/20*(-18*y3^4-48*y1*y3^3+80*y3^3*y2-80*y3^2*y2^2+5*h^2*y3^2+ ...
            120*y3^2*y2*y1-48*y3^2*y1^2-30*y3*y2*h^2+20*h^2*y3*y1- ...
            80*y3*y1*y2^2+120*y3*y2*y1^2-48*y3*y1^3+60*y2^2*h^2-78*y1^4+ ...
            35*h^2*y1^2-90*y2*h^2*y1+160*y2*y1^3-80*y1^2*y2^2)*V0/(y3-y1)/h^3;
        F2 = -1/10*(-28*y3^4-68*y1*y3^3+120*y3^3*y2-80*y3^2*y2^2+5*h^2*y3^2+ ...
            120*y3^2*y2*y1-48*y3^2*y1^2-68*y3*y1^3+50*h^2*y3*y1-60*y3*y2*h^2- ...
            80*y3*y1*y2^2+120*y3*y2*y1^2+5*h^2*y1^2-28*y1^4+60*y2^2*h^2- ...
            60*y2*h^2*y1-80*y1^2*y2^2+120*y2*y1^3)*V0/(y3-y1)/h^3;
        F3 = 1/20*(-78*y3^4-48*y1*y3^3+160*y3^3*y2-80*y3^2*y2^2+35*h^2*y3^2+ ...
            120*y3^2*y2*y1-48*y3^2*y1^2-90*y3*y2*h^2+20*h^2*y3*y1- ...
            80*y3*y1*y2^2+120*y3*y2*y1^2-48*y3*y1^3+60*y2^2*h^2-18*y1^4+ ...
            5*h^2*y1^2-30*y2*h^2*y1+80*y2*y1^3-80*y1^2*y2^2)*V0/(y3-y1)/h^3;
        % forces contains the equivalent nodal loads for
        % parabolic shear stress distribution at beam tip
        forces(2*(i-1)+1:2*i+1)=forces(2*(i-1)+1:2*i+1) + [F1 F2 F3]';
        % forcesb contains the equivalant nodal loads for
        % linear bending stress applied at the beam tip
        forcesb(2*(i-1)+1:2*i+1)=forcesb(2*(i-1)+1:2*i+1) + [FA FB FC]';
    end
end

if load_opt==3
    Interior = find(sqrt(coord(:,2).^2+(coord(:,3)-R).^2) - (R-h/2) < 0.001);
    LoadMat = [Interior(1) 2 -Mag/6];
    for i=2:(length(Interior)-1)
        dx = coord(Interior(i),2);
        dy = coord(Interior(i),3)-R;
        dL = sqrt(dx^2 + dy^2);
        if mod(i,2)==0
            Loadx = 2/3*Mag*dx/dL;
            Loady = 2/3*Mag*dy/dL;
        else
            Loadx = 1/3*Mag*dx/dL;
            Loady = 1/3*Mag*dy/dL;
        end
        LoadMat = [LoadMat
                   Interior(i) 1 Loadx];
        LoadMat = [LoadMat
                   Interior(i) 2 Loady];
    end
    LoadMat = [LoadMat
               Interior(end) 1 1/6*Mag];
end

if load_opt==4
    Interior = find(sqrt(coord(:,2).^2+(coord(:,3)-R).^2) - (R-h/2) < 0.001);
    Exterior = find(sqrt(coord(:,2).^2+(coord(:,3)-R).^2) > (R+h/2)-0.001);
    
    DispMat = [Interior(1) 2 -RadD];
    for i=2:(length(Interior)-1)
        dx = coord(Interior(i),2);
        dy = coord(Interior(i),3)-R;
        dL = sqrt(dx^2 + dy^2);
        Dx = RadD*dx/dL;
        Dy = RadD*dy/dL;
        DispMat = [DispMat
                   Interior(i) 1 Dx
                   Interior(i) 2 Dy];
    end
    DispMat = [DispMat
               Interior(end) 1 RadD];
end

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
