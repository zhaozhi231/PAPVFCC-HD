function [ dd ] = DIST(XS, X)

         
          dd = acosh(max(-LORENTZ(X, XS),1)) ;


end