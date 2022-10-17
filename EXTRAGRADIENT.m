function  [ X, iterk, f_eval, beta0, beta, rdist0, rdist, TimeCost ] = EXTRAGRADIENT(X, ITmax, Tolerence, rho)

% This is a Riemannian extragradient type method used to
% find the zeros of (X_1X_3, X_2X_3, X_3^2-1) on H^2.

                   delta = 1e-5 ; 
 
                   XS = [ 0 0 1]' ; 
            
%% Begining;

                    t0 =clock ;           
                    
          %  F(X) = (-X_1X3, -X_2X_3, 1-X_3^2) ;    
               
                    FX = FF(X) ;
                   
          %    The objective function  ;               
           
               beta0  = sqrt(LORENTZ(FX, FX)) ;

                    beta = beta0 ;    

                 rdist0 = DIST(XS, X) ;
     

%% Extragradient Method;

                   
                     XK = -FX ;

                       k = 0 ;   
                   
                     f_eval = 1 ;            
                     
% Stopping Criterion ;
                                        
              while   beta > Tolerence && k < ITmax                                             
                  
                        ak = 1 ;
   
                        
%% Intial Guess of The Step-size; 


                [ Z, FZ, lfro, rfro ] = PHI(X, XK, ak) ;        
             
                      f_eval = f_eval + 1 ;                                                      

                         l = 0 ;
                                       
                              while  lfro < delta*rfro ; 
                                                               
                                         ak = rho*ak ;
                                          
                                 [ Z, FZ, lfro, rfro ] = PHI(X, XK, ak) ;
                         
                                       f_eval = f_eval + 1 ;
                                       
                                       l = l + 1 ;
                                                                    
                              end
                              
%% Projection Step ;                         
                            
                       beta2 = LORENTZ(FZ,FZ)^(1/2) ; 
 
             if   beta2  > Tolerence      
                              
                              
                    FZ = FZ/sqrt(LORENTZ(FZ,FZ)) ;
            
                      aa = LORENTZ(FZ,FZ) ;
                      
                      bb = LORENTZ(X,FZ) ;
                       
                      PX = aa*X - bb*FZ ;
                      
                      if PX(3) > 0 ;
                          
                          mm = 1/sqrt(aa*(1+bb^2)) ;                          
                          
                      else
                          
                         mm = -1/sqrt(aa*(1+bb^2)) ;  

                      end
                      
                      
                          X = mm*PX ;
                         
                          FX = FF(X) ;  
           
                          XK = -FX ;
                          
          %    The objective function  ;               
           
                 beta  =  sqrt(LORENTZ(FX, FX)) ;
                 
%          fprintf('RSCG: Norm of (X_1X_3, X_2X_3, X_3^2-1) %d \n',beta)   
                                                       
             else
                 
                          beta = beta2 ;
                          
%          fprintf('RSCG: Norm of (X_1X_3, X_2X_3, X_3^2-1) %d \n',beta)  
                    
                          
                              X = Z ;
                          
                             FX = FZ ;
                 
                 
             end

                            k = k + 1 ;               
 
              end              
              
                       rdist = DIST(XS, X)  ; 
              
                            iterk = k ;    
  
                     TimeCost = etime(clock,t0) ;

                            
          fprintf('\n');
          fprintf('ETRGRAD: Number of Iterations %d \n', iterk)
          fprintf('ETRGRAD: Number of function estimations %d \n',f_eval)  
          fprintf('ETRGRAD: Initial norm of (X_1X_3, X_2X_3, X_3^2-1) ==========%d \n',beta0)
          fprintf('ETRGRAD: Final norm of (X_1X_3, X_2X_3, X_3^2-1) ==========%d \n',beta)
          fprintf('ETRGRAD: Initial Riemannian distance of X_k and X_* ==========%d \n',rdist0)
          fprintf('ETRGRAD: Final Riemannian distance of X_k and X_* ==========%d \n',rdist)
          fprintf('ETRGRAD: Computing time used ========== %d \n', TimeCost)                               
         

end

