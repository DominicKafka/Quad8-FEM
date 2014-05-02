function [Res,Tangent,stress,strain] = Quad8_Res_and_Tangent(X,U,Cmat,t)
Tangent   = zeros(16,16);        % Initialize tangent with zeros
Res       = zeros(16,1);         % Initialize residual with zeros
Gauss_pos = 1/sqrt(3);           % Gauss point location
Ivec      = [1, 1, 0, 0]';       % Identity vector Eq.(4.89)
stress    = zeros(12,1);         % Initialize stress vector with zeros
strain    = zeros(12,1);         % Initialize strain vector with zeros
for jGauss = 1:2                 % 2 by 2 Gauss integration loops
    eta = (-1)^jGauss*Gauss_pos; % Natural coordinate eta
    for iGauss = 1:2
        xi = (-1)^iGauss*Gauss_pos;      
		% Natural coordinate xi
        [B,detJ] = B_Quad8(xi,eta,X,true);      % B for Eq.(4.83)
        Fvec = Ivec + B*U;                      % Eq.(4.88)
        %detF = Fvec(1)*Fvec(2)-Fvec(3)*Fvec(4); % Determinant of F
        Fmat = Fvec_to_Fmat(Fvec);              % Eq.(4.86)
% Mooney-Rivlin case
%        Cvec = Fmat'*Fvec;
%        [Svec,Cmat] = MooneyRivlin(Cvec);
% Saint Venant- Kirchhoff case
        Evec = 0.5*(Fmat'*Fvec - Ivec);         % Eq.(4.90)
        Svec = Cmat*Evec;                       % Eq.(4.96)
        Smat = Svec_to_Smat(Svec);              % Eq.(4.87)
        Res  = Res + t*B'*Fmat*Svec*detJ;       % Eq.(4.98)
        Tangent = Tangent + t*B'*(Smat + Fmat*Cmat*Fmat')*B*detJ; % Eq.(4.106)
        detF = Fvec(1)*Fvec(2)-Fvec(3)*Fvec(4);
        FTF = [Fvec(1)^2      Fvec(3)^2    Fvec(1)*Fvec(3) Fvec(3)*Fvec(1)
               Fvec(4)^2      Fvec(2)^2    Fvec(4)*Fvec(2) Fvec(4)*Fvec(2)
               Fvec(1)*Fvec(4) Fvec(3)*Fvec(2) Fvec(1)*Fvec(2) Fvec(3)*Fvec(4)];
        Cauchy                   = FTF*Svec/detF;
        GP                       = 2*(jGauss-1) + iGauss;
        stress(3*(GP-1)+(1:3),1) = Cauchy;
        strain(3*(GP-1)+(1:3),1) = Evec(1:3);
    end
end
