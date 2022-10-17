function [ Z, FZ, lfro, rfro ] = PHI2(X, XK, PCBF, ak, bk)
% This is the algorithm used to calculate
%    F(X) = (-X_1X3, -X_2X_3, 1-X_3^2) ;


         cc =  sqrt(LORENTZ(XK,XK)) ;

             XK = XK/cc ;
             
             ak = cc*ak ;
        
        VX = sinh(ak)*X + cosh(ak)*XK ;
        
        Z = cosh(ak)*X + sinh(ak)*XK ;    

               FZ = FF(Z)  ;


          lfro = -LORENTZ(FZ, VX) ;
  
         rfro = ak/bk*(DIST(X,PCBF)^2) ;
        

end

