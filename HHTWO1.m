function  HHTWO1

% This is the comparision of different kinds of numerical methods for
% solving mnonotone vector field on Hadamard manifolds.
   

%% Parameter Settings ;
  
             ITmax = 1e3;     Tolerence = 1e-6 ;
    
                        rho = 0.5 ;
                     
      % Rpeated number of numerical experiments for the same problem size;

                          T = 10 ;

             

% Extragradient Method;

% Average number of iteration;
IT1  = 0 ;
% Average number of function evaluation ;
NF1 = 0 ;
% Average value of initial norm of (X_1X_3, X_2X_3, X_3^2-1) ;
INMFX1 = 0 ;
% Average value of final norm of (X_1X_3, X_2X_3, X_3^2-1) ;
FNMFX1 = 0 ;
% Average value of initial Riemannian distance of X_0 and X_* ;
IRIDST1 = 0 ;
% Average value of final Riemannian distance of X_k and X_* ;
FRIDST1 = 0 ;
% Average cputime ;
TIME1 = 0 ;



% Riemannian INERTIAL Method;

% Average number of iteration;
IT2  = 0 ;
% Average number of function evaluation ;
NF2 = 0 ;
% Average value of initial norm of (X_1X_3, X_2X_3, X_3^2-1) ;
INMFX2 = 0 ;
% Average value of final norm of (X_1X_3, X_2X_3, X_3^2-1) ;
FNMFX2 = 0 ;
% Average value of initial Riemannian distance of X_0 and X_* ;
IRIDST2 = 0 ;
% Average value of final Riemannian distance of X_k and X_* ;
FRIDST2 = 0 ;
% Average cputime ;
TIME2 = 0 ;


l = 1 ;

while l <= T   
                   
                               
%% Initial Guesses; 
                 
                   X = rand(2,1) ;
                
               X = [ X ; sqrt(X'*X+1) ] ;     
                 
                 
                   
%% Begining;
                         
% Extragradient Method;

  [ X1, iterk1, f_eval1, Inbeta1, beta1, rdist01, rdist1, TimeCost1 ] = EXTRAGRADIENT(X, ITmax, Tolerence, rho) ; 
  

IT1  = IT1 + iterk1 ;
NF1 = NF1 + f_eval1 ;
INMFX1 = INMFX1 + Inbeta1 ;
FNMFX1 = FNMFX1 + beta1 ;
IRIDST1 = IRIDST1 + rdist01 ;
FRIDST1 = FRIDST1 + rdist1 ;
TIME1 = TIME1 + TimeCost1 ;

 
% Riemannian INERTIAL Method;


    [ X5, iterk2, f_eval2, Inbeta2, beta2, rdist02, rdist2, TimeCost2 ] = INERTIAL1(X, ITmax, Tolerence, rho) ;  
            
             
IT2  = IT2 + iterk2 ;
NF2 = NF2 + f_eval2 ;
INMFX2 = INMFX2 + Inbeta2 ;
FNMFX2 = FNMFX2 + beta2 ;
IRIDST2 = IRIDST2 + rdist02 ;
FRIDST2 = FRIDST2 + rdist2 ;
TIME2 = TIME2 + TimeCost2 ;

   
l = l + 1 ;

end


IT1 = IT1/T; 
NF1 = NF1/T; 
INMFX1 = INMFX1/T;
FNMFX1 = FNMFX1/T;
IRIDST1 = IRIDST1/T;
FRIDST1 = FRIDST1/T;
TIME1 = TIME1/T;


IT2 = IT2/T; 
NF2 = NF2/T;
INMFX2 = INMFX2/T;
FNMFX2 = FNMFX2/T;
IRIDST2 = IRIDST2/T;
FRIDST2 = FRIDST2/T;
TIME2 = TIME2/T;


TIME1
TIME2

IT1
IT2


NF1
NF2


INMFX1
INMFX2

FNMFX1
FNMFX2


IRIDST1
IRIDST2


FRIDST1
FRIDST2




end

