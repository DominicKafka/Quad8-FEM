function [StressNode,VonMises,Tresca] = write_output_file(file_out,U,displ,Pb,plane,pois,nodes,elnodes,stress,strain)

ndispl = size(displ, 1);
nnodes = length(nodes);
nelem = size(elnodes, 1);

Uoutput = zeros(nnodes,3);
Uoutput(:,1) = nodes;
Uoutput(:,2) = U(1:2:nnodes*2);
Uoutput(:,3) = U(2:2:nnodes*2);

[StressOut, StressNode] = nodal_stresses(elnodes, stress);
VonMises = calc_von_mises(StressNode, pois, plane);
Tresca = calc_tresca(StressNode, pois, plane);

StrainOut = zeros(nelem,13);
StrainOut(1:nelem,1) = elnodes(:,1);
StrainOut(:,2:7) = strain(:,1:6);
StrainOut(:,8:10) = strain(:,10:12);
StrainOut(:,11:13) = strain(:,7:9);

SupReac = zeros(ndispl,3);
SupReac(:,1:2) = displ(:,1:2);
SupReac(:,3) = Pb;

fid = fopen(file_out,'w');
fprintf(fid,'OUTPUT OF MATLAB Q4 SMALL STRAIN FEM IMPLEMENTATION \n');
fprintf(fid,'\n');
fprintf(fid,'          DISPLACEMENTS \n');
fprintf(fid,'********************************* \n');
fprintf(fid,'  Node      U1           U2 \n');
fprintf(fid,'********************************* \n');
fprintf(fid,'%5d %13.5e %13.5e \n',Uoutput');

fprintf(fid,'\n');
fprintf(fid,'                ELEMENT STRESSES \n');
fprintf(fid,['*********************************************************', ...
             '*********************************************************', ...
             '*********************************************** \n']);
fprintf(fid,['Element   S11_G1       S22_G1       S12_G1       S11_G2   ', ...
             '    S22_G2       S12_G2       S11_G3       S22_G3       ', ...
             '  S12_G3     S11_G4       S22_G4       S12_G4 \n']);
fprintf(fid,['*********************************************************', ...
             '*********************************************************', ...
             '*********************************************** \n']);
fprintf(fid,['%5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e ', ...
             '%12.4e %12.4e %12.4e %12.4e %12.4e \n'],StressOut');

fprintf(fid,'\n');
fprintf(fid,'                ELEMENT STRAINS \n');
fprintf(fid,['*********************************************************', ...
             '*********************************************************', ...
             '*********************************************** \n']);
fprintf(fid,['Element   E11_G1       E22_G1       E12_G1       E11_G2   ', ...
             '    E22_G2       E12_G2       E11_G3       E22_G3     ', ...
             '  E12_G3       E11_G4       E22_G4       E12_G4 \n']);
fprintf(fid,['*********************************************************', ...
             '*********************************************************', ...
             '*********************************************** \n']);
fprintf(fid,['%5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e ', ...
             '%12.4e %12.4e %12.4e %12.4e %12.4e \n'],StrainOut');

fprintf(fid,'\n');
fprintf(fid,'       SUPPORT REACTIONS \n');
fprintf(fid,'***************************** \n');
fprintf(fid,'  Node   Dof      Magnitude \n');
fprintf(fid,'***************************** \n');
fprintf(fid,'%5d %5d %17.5e \n',SupReac');
fclose(fid);

finish = toc;
disp(['Done writing output                 : ',num2str(finish),' seconds.'])