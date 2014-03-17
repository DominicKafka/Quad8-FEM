function [StressNode,VonMises,Tresca] = write_output_file(file_out,U,displ,ndispl,Pb,plane,pois,elas,nnodes,nodes,nelem,elnodes,stress,strain)

Uoutput = zeros(nnodes,3);
Uoutput(:,1) = nodes;
Uoutput(:,2) = U(1:2:nnodes*2);
Uoutput(:,3) = U(2:2:nnodes*2);

StressOut = zeros(nelem,13);
StressOut(:,1)     = elnodes(:,1);
StressOut(:,2:7)   = stress(:,1:6);
StressOut(:,8:10)  = stress(:,10:12);
StressOut(:,11:13) = stress(:,7:9);

StressNode(1:nelem,1) = elnodes(:,1);
StressNode(:,2) = (1+sqrt(3)/2)*StressOut(:,2) - 0.5*(StressOut(:,5) + ...
                  StressOut(:,11)) + (1-sqrt(3)/2)*StressOut(:,8);
StressNode(:,5) = (1+sqrt(3)/2)*StressOut(:,5) - 0.5*(StressOut(:,2) + ...
                  StressOut(:,8)) + (1-sqrt(3)/2)*StressOut(:,11);
StressNode(:,8) = (1+sqrt(3)/2)*StressOut(:,8) - 0.5*(StressOut(:,5) + ...
                  StressOut(:,11)) + (1-sqrt(3)/2)*StressOut(:,2);
StressNode(:,11) = (1+sqrt(3)/2)*StressOut(:,11) - 0.5*(StressOut(:,2) + ...
                  StressOut(:,8)) + (1-sqrt(3)/2)*StressOut(:,5);

StressNode(:,3) = (1+sqrt(3)/2)*StressOut(:,3) - 0.5*(StressOut(:,6) + ...
                  StressOut(:,12)) + (1-sqrt(3)/2)*StressOut(:,9);
StressNode(:,6) = (1+sqrt(3)/2)*StressOut(:,6) - 0.5*(StressOut(:,3) + ...
                  StressOut(:,9)) + (1-sqrt(3)/2)*StressOut(:,12);
StressNode(:,9) = (1+sqrt(3)/2)*StressOut(:,9) - 0.5*(StressOut(:,6) + ...
                  StressOut(:,12)) + (1-sqrt(3)/2)*StressOut(:,3);
StressNode(:,12) = (1+sqrt(3)/2)*StressOut(:,12) - 0.5*(StressOut(:,3) + ...
                  StressOut(:,9)) + (1-sqrt(3)/2)*StressOut(:,6);

StressNode(:,4) = (1+sqrt(3)/2)*StressOut(:,4) - 0.5*(StressOut(:,7) + ...
                  StressOut(:,13)) + (1-sqrt(3)/2)*StressOut(:,10);
StressNode(:,7) = (1+sqrt(3)/2)*StressOut(:,7) - 0.5*(StressOut(:,4) + ...
                  StressOut(:,10)) + (1-sqrt(3)/2)*StressOut(:,13);
StressNode(:,10) = (1+sqrt(3)/2)*StressOut(:,10) - 0.5*(StressOut(:,7) + ...
                  StressOut(:,13)) + (1-sqrt(3)/2)*StressOut(:,4);
StressNode(:,13) = (1+sqrt(3)/2)*StressOut(:,13) - 0.5*(StressOut(:,4) + ...
                  StressOut(:,10)) + (1-sqrt(3)/2)*StressOut(:,7);

VonMises = zeros(nelem,4);
if (plane==1)
    for i=1:4
        VonMises(:,i) = StressNode(:,2+(i-1)*3).^2 - ...
                StressNode(:,2+(i-1)*3).*StressNode(:,3+(i-1)*3) + ...
                StressNode(:,3+(i-1)*3).^2 + 3*StressNode(:,4+(i-1)*3).^2;
        VonMises(:,i) = VonMises(:,i).^0.5;
    end
else
    for i=1:4
        VonMises(:,i) = (1-pois+pois^2)*(StressNode(:,2+(i-1)*3).^2 + ...
                        StressNode(:,3+(i-1)*3).^2) - (1+pois-pois^2)* ...
                        StressNode(:,2+(i-1)*3).*StressNode(:,3+(i-1)*3)+...
                        3*StressNode(:,4+(i-1)*3).^2;
        VonMises(:,i) = VonMises(:,i).^0.5;
    end
end    

Tresca = zeros(nelem,4);
if plane ==1
    for j=1:nelem
        for i=1:4
            s = [StressNode(j,2+(i-1)*3) StressNode(j,4+(i-1)*3) 0
                 StressNode(j,4+(i-1)*3) StressNode(j,3+(i-1)*3) 0
                0 0 0];
            principal = eig(s);
            Tresca(j,i) = max(principal)-min(principal);
        end
    end
else
    for j=1:nelem
        for i=1:4
            s = [StressNode(j,2+(i-1)*3) StressNode(j,4+(i-1)*3) 0
                 StressNode(j,4+(i-1)*3) StressNode(j,3+(i-1)*3) 0
                0 0 pois*(StressNode(j,2+(i-1)*3)+StressNode(j,3+(i-1)*3))];
            principal = eig(s);
            Tresca(j,i) = max(principal)-min(principal);
        end
    end
end
       
StrainOut = zeros(nelem,13);
StrainOut(1:nelem,1) = elnodes(:,1);
StrainOut(:,2:7) = strain(:,1:6);
StrainOut(:,8:10) = strain(:,10:12);
StrainOut(:,11:13) = strain(:,7:9);

SupReac = zeros(ndispl,3);
SupReac(:,1:2) = displ(:,1:2);
SupReac(:,3) = Pb;

fid = fopen([file_out],'w');
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