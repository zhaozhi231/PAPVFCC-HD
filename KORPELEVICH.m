function [ X, iterk, f_eval, beta0, beta, rdist0, rdist, TimeCost ] = KORPELEVICH(X, ITmax, Tolerence, rho)

% This is a Riemannian Korpelevich type method used to
% find the solution of a variational ineuality defined
% on a convex subset of H^2.

%% Parameter Settings ;
      
            delta = 1e-5 ;   bk = 1 ;                      
                                                   
                 XS = [ 0 0 1]' ;   
                 
%% Begining;

                    t0 =clock ;           
                    
          %  F(X) = (-X_1X3, -X_2X_3, 1-X_3^2) ;    
               
                    FX = FF(X) ;
                    
                 BFX = EXP(X,-FX,bk) ;
                 
                  PCBF = PC(BFX) ;
                 
                RXB = EXPINV(X,PCBF) ;
                   
          %    The objective function  ;               
           
             beta0  = sqrt(LORENTZ(FX,FX)) ;

                     beta = beta0 ;    
                     
                 rdist0 = DIST(XS, X) ;   
                   
                       XK = RXB ;
     

%% Extragradient Method;
          
                       k = 0 ;   
                   
                     f_eval = 1 ;            
                     
% Stopping Criterion ;
                                        
              while   beta > Tolerence && k < ITmax                                             
                  
                        ak = 1 ;
   
                        
%% Intial Guess of The Step-size; 


                [ Z, FZ, lfro, rfro ] = PHI2(X, XK, PCBF, ak, bk) ;           
             
                      f_eval = f_eval + 1 ;                                                      
    
                              while  lfro < delta*rfro ; 
                                                               
                                         ak = rho*ak ;
                                          
                                 [ Z, FZ, lfro, rfro ] = PHI2(X, XK, PCBF, ak, bk) ;
                         
                                       f_eval = f_eval + 1 ;
                                                        
                              end
                              
%% Projection Step ;          
                              
                              
                    FZ = FZ/sqrt(LORENTZ(FZ,FZ)) ;
            
                      aa = LORENTZ(FZ,FZ) ;
                      
                      bb = LORENTZ(X,FZ) ;
                       
                      PX = aa*X - bb*FZ ;
                      
                      if PX(3) > 0 ;
                          
                          mm = 1/sqrt(aa*(1+bb^2)) ;                          
                          
                      else
                          
                         mm = -1/sqrt(aa*(1+bb^2)) ;  

                      end
                      
                      
                           W = mm*PX ;         
                         
                           X = PC(W) ;
                         
                           FX = FF(X) ;
                    
                       BFX = EXP(X,-FX,bk) ;
                 
                        PCBF = PC(BFX) ;
                 
                      RXB = EXPINV(X,PCBF) ;
                   
               %    The objective function  ;               
           
                  beta  = sqrt(LORENTZ(FX,FX)) ;

                          XK = RXB ;
                                          
                           k = k + 1 ;
                                 

              end              
              
                       rdist = DIST(XS, X) ;   
              
                            iterk = k ;    
  
                     TimeCost = etime(clock,t0) ;

                            
          fprintf('\n');
          fprintf('KORPELEVICH: Number of Iterations %d \n', iterk)
          fprintf('KORPELEVICH: Number of function estimations %d \n',f_eval)  
          fprintf('KORPELEVICH: Initial norm of (X_1X_3, X_2X_3, X_3^2-1) ==========%d \n',beta0)
          fprintf('KORPELEVICH: Final norm of (X_1X_3, X_2X_3, X_3^2-1) ==========%d \n',beta)
          fprintf('KORPELEVICH: Riemannian distance of X_0 and X_* ==========%d \n',rdist)
          fprintf('KORPELEVICH: Riemannian distance of X_k and X_* ==========%d \n',rdist)
          fprintf('KORPELEVICH: Computing time used ========== %d \n', TimeCost)                
                          
         

end

