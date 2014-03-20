clear all;close all
%ElType = '5B';
%ElType = '7B';
%ElType = 'q4';
ElType = 'q8';

if any(strcmp(ElType, {'q4', '5B', '7B'}))
    NodesPerEl = 4;
elseif strcmp(ElType, 'q8')
    NodesPerEl = 8;
end
DofPerEl   = 2*NodesPerEl;
TotKvec    = DofPerEl^2;

disp('Welcome to the MEE781 Finite Element Program')

filename = input('Input filename (without extension) ? ','s');
[nnodes,ndcoor,nodes,coor,nelem,plane,elnodes,elas,pois,t,ndispl,...
    displ,ncload,cload,nloadinc,mdof,sdof] = read_input_file(filename,NodesPerEl);
file_out = [filename,'.out'];

%Get maximum dimensions of model
dx_max = max(coor(:,1))-min(coor(:,1));
dy_max = max(coor(:,2))-min(coor(:,2));
dL_max = sqrt(dx_max^2 + dy_max^2);

%GraphOpt=1: graphical display of results and text based output files
%GraphOpt=0: only text based output files
GraphOpt = false;

% Find prescribed (pdof) and free (fdof) degrees of freedom
dof = ones(nnodes*2,1);
Up  = zeros(ndispl,1);
for i=1:ndispl
    pos=(find(nodes==displ(i,1))-1)*2+displ(i,2);
    dof(pos,1)=0;
    Up(i,1)=displ(i,3);
end
pdof=find(dof==0);
fdof=find(dof~=0);
fdof=setdiff(fdof,mdof);
fdof=setdiff(fdof,sdof);

% Initially guess that all free displacements are zero
U = zeros(2*nnodes,1);
%U(fdof,1) = zeros(length(fdof),1);

% Construct elasticity tensor C
% If plane strain
if plane==0
    e = elas/(1-pois^2);
    nu = pois/(1-pois);
else
    e = elas;
    nu = pois;
end

c = e/(1-nu^2);
matC = zeros(4,4);
matC(1,1) = c;
matC(2,2) = c;
matC(1,2) = c*nu;
matC(2,1) = matC(1,2);
matC(3,3) = c*(1-nu);
matC(4,4) = matC(3,3);
ndnum  = 2:(1+NodesPerEl);
rowpos = repmat(1:(DofPerEl),1,DofPerEl);
colpos = reshape(repmat(1:DofPerEl,DofPerEl,1),1,TotKvec);

tol    = 3.d-5;
dUNrm  = 1.0;

F = zeros(2*nnodes,1);
LoadFac  = (1:1:nloadinc)/nloadinc;

% FIXME: this 12 should be a variable
stress = zeros(nelem, 12);
strain = zeros(nelem, 12);

for iter_load = 1:nloadinc;
    disp('---------------------------------------------------------------')
    disp(['                      Load increment ',num2str(iter_load)])
    disp('---------------------------------------------------------------')
    % Correct fraction of prescribed displacement applied to pdof
    U(pdof,1) = LoadFac(iter_load)*Up;
    % If MPCs are present, compute slave DOF values
    if ~isempty(sdof)
        [U(sdof,1),P_mpc,d2PdUm2_Rs] = MPC_user(U(mdof,1),F(sdof));
    else
        P_mpc      = [];
        d2PdUm2_Rs = [];
    end
    ResNrm = 1.d0;
    iter   = 0;
    while (ResNrm>tol) || (dUNrm>tol)
        tic;
        iter = iter + 1;
        %Main loop over elements. Compute k_elem and assemble
        %Initialize global stiffness matrix vectors;
        row_vec   = zeros(TotKvec*nelem,1);
        col_vec   = zeros(TotKvec*nelem,1);
        stiff_vec = zeros(TotKvec*nelem,1);
        Residual  = zeros(2*nnodes,1);
        %Initialize global load vector
        F_ext = zeros(2*nnodes,1);
        pos_vec = 0;
        for i=1:nelem
            % Find reference coordinates of element nodes
            X = coor(elnodes(i,ndnum),1);
            Y = coor(elnodes(i,ndnum),2);
            % Get global degree of freedom numbers per element
            pg          = zeros(DofPerEl,1);
            pg(1:2:end) = 2*elnodes(i,2:(1+NodesPerEl))-1;
            pg(2:2:end) = 2*elnodes(i,2:(1+NodesPerEl));
            % Get current guess for nodal displacements
            U_el = U(pg);
            XY   = [X Y];
            switch ElType
                case 'q4'
                    [El_res,k_elem,El_stress,El_strain] = Quad4_Res_and_Tangent(XY,U_el,matC,t);
                case '5B'
                    [El_res,k_elem,El_stress,El_strain] = FiveB_Res_and_Tangent(XY,U_el,matC,t);
                case '7B'
                    [El_res,k_elem,El_stress,El_strain] = SevenB_Res_and_Tangent(XY,U_el,matC,t);
                case 'q8'
                    [El_res,k_elem,El_stress,El_strain] = Quad8_Res_and_Tangent(XY,U_el,matC,t);
            end
            stress(i,:) = El_stress;
            strain(i,:) = El_strain;
            % Assemble residual
            Residual(pg) = Residual(pg) + El_res;
            % Assemble k_elem into sparse k_global using vectors
            row_vec(TotKvec*(i-1)+(1:TotKvec),1)   = pg(rowpos);
            col_vec(TotKvec*(i-1)+(1:TotKvec),1)   = pg(colpos);
            stiff_vec(TotKvec*(i-1)+(1:TotKvec),1) = k_elem(1:TotKvec);
        end;  % End of main loop over elements
        % Assemble k_global from vectors
        k_global = sparse(row_vec,col_vec,stiff_vec,2*nnodes,2*nnodes);
        clear row_vec col_vec stiff_vec
        finish = toc;
        disp(['Done assembling stiffness matrix: ',num2str(finish),' seconds.'])
        
        % Add nodal loads to global load vector
        for i=1:ncload;
            p = find(nodes==cload(i,1));
            pos = (p-1)*2+cload(i,2);
            F_ext(pos,1) = F_ext(pos,1)+LoadFac(iter_load)*cload(i,3);
        end        
        
        % Subtract internal nodal loads
        F = Residual - F_ext;
        
        ResNrm = norm(F(fdof,1));
        if iter == 1
            if ResNrm > 1e-4
                ResNrm0 = ResNrm;
            else
                ResNrm0 = 1;
            end
        end
        ResNrm = ResNrm/ResNrm0;
        PrntResStr = ['Normalized residual at start of iteration ',num2str(iter),'    = ',num2str(ResNrm,'%10.6e')];
        disp(PrntResStr)
        
        % Solve update of free dof's
        % Solution for non-symmetric stiffness matrix
        
        Kff = k_global(fdof,fdof);
        Pf  = F(fdof);
        Kfm = k_global(fdof,mdof) + k_global(fdof,sdof)*P_mpc;
        Kmm = k_global(mdof,mdof) + k_global(mdof,sdof)*P_mpc + ...
              P_mpc'*k_global(sdof,mdof) + P_mpc'*k_global(sdof,sdof)*P_mpc + ...
              d2PdUm2_Rs;
        % Define RHS
        Pm  = F(mdof) + P_mpc'*(F(sdof));
        Pa  = [Pf
               Pm];
        Kaa = [Kff  Kfm
               Kfm' Kmm];
        
        finish = toc;
        disp(['Done assembling stiffness matrix: ',num2str(finish),' seconds.'])
        
        tic;
        
        deltaUf = -Kaa\Pa;
        
        finish  = toc;
        disp(['Done solving system             : ',num2str(finish),' seconds.'])
        
        dUNrm = norm(deltaUf)/dL_max;
        PrntDspStr = ['Normalized displacement update                 = ',num2str(dUNrm,'%10.6e')];
        disp(PrntDspStr)
        disp('                    --------------------')
        % Sort Uf and Ub into A
        U(fdof,1) = U(fdof,1) + deltaUf(1:length(fdof));
        if ~isempty(sdof)
            U(mdof,1) = U(mdof,1) + deltaUf(length(fdof)+1:end);
            [U(sdof,1),P_mpc,d2PdUm2_Rs] = MPC_user(U(mdof,1),F(sdof));
        end
        AllResNrm(iter) = ResNrm; %#ok<SAGROW>
        AlldUNrm(iter)  = dUNrm; %#ok<SAGROW>
    end
    % Get support reactions
    Fp = F(pdof);

    disp(['Load increment ',num2str(iter_load),' converged after ',num2str(iter),' iterations.'])
    All_iter(iter_load) = iter; %#ok<SAGROW>
    All_soln(1+nnodes*(iter_load-1):nnodes*iter_load,:) = [U(1:2:2*nnodes) U(2:2:2*nnodes)];
end
tic;


% Compute Von Mises and Tresca and write output to text based output file
[StressNode,VonMises,Tresca] = write_output_file(file_out,U,displ,Fp, ...
    plane,pois,nodes,elnodes,stress,strain);

% If GraphOpt=1, start Graphical Output
if GraphOpt
    graphical_user_interface(nnodes,coor,nelem,elnodes,StressNode,U,VonMises,Tresca,nloadinc,All_soln);
end


