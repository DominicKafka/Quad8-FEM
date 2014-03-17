function [B,detJ] = B_Quad8(xi,eta,X,NL_flag)
D     = zeros(4,16); % Initialize D with zeros
% Derivatives of shape functions wrt xi & eta
dNdxi = [1/4*eta+1/2*xi*(1-eta)-1/4*eta^2,-1/4*eta+1/2*xi*(1-eta)+1/4*eta^2,...
         1/4*eta+1/4*eta^2+1/2*xi*(1+eta),-1/4*eta+1/2*xi*(1+eta)-1/4*eta^2,...
         -xi*(1-eta),1/2-1/2*eta^2, -xi*(1+eta),-1/2+1/2*eta^2
         1/4*xi-1/4*xi^2+(1/2-1/2*xi)*eta,-1/4*xi-1/4*xi^2+(1/2+1/2*xi)*eta,...
         1/4*xi+(1/2+1/2*xi)*eta+1/4*xi^2,-1/4*xi+1/4*xi^2+(1/2-1/2*xi)*eta,...
         -1/2+1/2*xi^2,-2*(1/2+1/2*xi)*eta,1/2-1/2*xi^2,-2*(1/2-1/2*xi)*eta];
D(1:2,1:2:15) = dNdxi;
D(3:4,2:2:16) = dNdxi; % Arrange shape function derivatives into D
J      = dNdxi*X;      % Eq.(2.40)
detJ   = det(J);       % Determinant of Jacobian J
invJ   = inv(J);       % Inverse of Jacobian
if NL_flag             % Four rows required for nonlinear case
    dxidX  = [invJ(1,1), invJ(1,2),    0     ,    0
                 0     ,    0     , invJ(2,1), invJ(2,2)
              invJ(2,1), invJ(2,2),     0     ,    0
                 0     ,    0     , invJ(1,1), invJ(1,2)];
else                  % Three rows required in linear case
    dxidX  = [invJ(1,1), invJ(1,2),    0     ,    0
                 0     ,    0     , invJ(2,1), invJ(2,2)
              invJ(2,1), invJ(2,2), invJ(1,1), invJ(1,2)];
end
B = dxidX*D; %Shape function derivatives wrt x and y: Eq.(2.39)
