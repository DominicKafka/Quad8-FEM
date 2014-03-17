function [Smat] = Svec_to_Smat(Svec);
Smat = [Svec(1) ,  0      , Svec(3) ,  0
         0      , Svec(2) ,  0      , Svec(3)
        Svec(3) ,  0      , Svec(2) ,  0
         0      , Svec(3) ,  0      , Svec(1)]; % Eq.(4.87)