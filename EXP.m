function [ NX ] = EXP(X, XK, ak)
%  This algorithm calculate the retraction of vectors
%  onto the tangent spaves;
%             R_(X)(ak*XK)

%% R_X(ak*XK);

     
          cc =  sqrt(LORENTZ(XK,XK)) ;

                  XK = XK/cc ;
                               
                  ak = cc*ak ;

        NX = cosh(ak)*X + sinh(ak)*XK ;
        
         NX = NX/sqrt(-LORENTZ(NX,NX)) ;
        

end

