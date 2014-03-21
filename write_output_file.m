function write_output_file(file_out,U,displ,Pb,nodes,elnodes,strain,StressOut)

Uoutput = [nodes, U(1:2:end), U(2:2:end)];

StrainOut = [elnodes(:,1), strain(:, [1:6, 10:12, 7:9])];

SupReac = [displ(:,1:2), Pb];

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