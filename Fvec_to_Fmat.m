function [Fmat] = Fvec_to_Fmat(Fvec);
Fmat = [Fvec(1),  0      , 0.5*Fvec(3) , 0.5*Fvec(3)
         0     , Fvec(2) , 0.5*Fvec(4) , 0.5*Fvec(4)
         0     , Fvec(3) , 0.5*Fvec(1) , 0.5*Fvec(1)
        Fvec(4),  0      , 0.5*Fvec(2) , 0.5*Fvec(2)]; % Eq.(4.86)