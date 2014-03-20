function [nnodes,ndcoor,nodes,coor,nelem,plane,elnodes,elas,pois,t,...
    ndispl,displ,ncload,cload,nloadinc,MasterDOF,SlaveDOF] = read_input_file(filename,NodesPerEl)

% Open input file
fid = fopen([filename,'.inp'],'r');

tic;
% Read number of nodes
dummy = fgetl(fid);
nnodes = fscanf(fid,'%d \n',1);

%Read nodal coordinates
dummy = fgetl(fid);
%dummy = fgetl(fid);
ndcoor = fscanf(fid,'%f \n',[3,nnodes])';

% Store node numbers in nodes
nodes = round(ndcoor(:,1));

% Store coordinates in coord
coor = ndcoor(:,2:3);

% Read number of elements
dummy = fgetl(fid);
nelem = fscanf(fid,'%5d \n',1);

% Read if elements are plane stress or plain strain
dummy = fgetl(fid);
plane = fscanf(fid,'%1d \n',1);

% Read element number and element connectivity
dummy = fgetl(fid);
elnodes = fscanf(fid,'%5d \n',[1+NodesPerEl,nelem])';

% Read material constants
dummy = fgetl(fid);
elas = fscanf(fid,'%f',1);
pois = fscanf(fid,'%f',1);

% Read element thickness
t = fscanf(fid,'%f \n',1);

% Read number of prescribed displacements
dummy = fgetl(fid);
ndispl = fscanf(fid,'%5d \n',1);

% Read prescribed displacements;
dummy = fgetl(fid);
displ = fscanf(fid,'%f \n',[3,ndispl])';
displ = sortrows(displ);

% Read number of nodal loads
dummy = fgetl(fid);
ncload = fscanf(fid,'%5d \n',1);

% Read nodal loads
dummy = fgetl(fid);
cload = fscanf(fid,'%f \n',[3,ncload])';

% Read number of load increments
dummy = fgetl(fid);
nloadinc = fscanf(fid,'%d \n',1);

% Read number of MPCs
dummy = fgetl(fid);
nMPC = fscanf(fid,'%5d \n',1);

if nMPC>0
    % Read number of master dofs in MPCs
    dummy = fgetl(fid);
    nMaster = fscanf(fid,'%5d \n',1);
    
    %Read master dofs
    dummy = fgetl(fid);
    MasterDOF = zeros(nMaster,2);
    for i=1:nMaster
        MasterDOF(i,:) = fscanf(fid,'%5d %5d \n',2);
    end
    
    % Read slave DOFs
    dummy = fgetl(fid);
    SlaveDOF = zeros(nMPC,2);
    for i=1:nMPC
        SlaveDOF(i,:) = fscanf(fid,'%5d %5d \n',2);
    end
    
    MasterDOF = 2*(MasterDOF(:,1)-1) + MasterDOF(:,2);
    SlaveDOF  = 2*(SlaveDOF(:,1)-1) + SlaveDOF(:,2);
    
else
    MasterDOF = [];
    SlaveDOF  = [];
end
fclose(fid);
finish = toc;
disp(['Done reading input file             : ',num2str(finish),' seconds'])
