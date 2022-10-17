function  [ X, iterk, f_eval, beta0, beta, rdist0, rdist, TimeCost ] = INERTIAL1(X, ITmax, Tolerence, rho)


% This is a Riemannian Inertial type method used to
% find the zeros of (X_1X_3, X_2X_3, X_3^2-1) on H^2.

%% Parameter Settings ; 

                   sigma = 1e-6 ;
            
               nu = 0.1 ;         mu = 1e-2 ;
               
                   XS = [ 0 0 1]' ; 
                 
%% Begining;

                    t0 =clock ;           
                    
          %  F(X) = (-X_1X3, -X_2X_3, 1-X_3^2) ;    
               
                    FX = FF(X) ;
                   
          %    The objective function  ;               
           
               beta0  = sqrt(LORENTZ(FX, FX)) ;

                    beta = beta0 ;
                    
      fprintf('INERTIAL: Norm of (X_1X_3, X_2X_3, X_3^2-1) %d \n',beta)  
                    
                 rdist0 = DIST(XS, X) ;
                    
                    
%% Initial Nonlinear Conjugate Direction;    
                    
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
                                       
                              while  lfro < sigma*rfro ; 
                                                               
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
                      

                            NX = mm*PX ;

                           NFX = FF(NX) ;   
                         
                         
%%  Inertial Strategy;



                           dd = LORENTZ(NFX,NFX) ;
                             
          %    The objective function  ;               
               
                             beta = sqrt(dd) ;
          
%             fprintf('INERTIAL: Norm of (X_1X_3, X_2X_3, X_3^2-1) %d \n',beta)   
                       
 
                             ss = EXPINV(X,NX) ;
                          
                           aa = LORENTZ(ss,FX) ;    
                           
                            nu = 0.1/(k+3) ;

                       bk = nu*dd/(abs(aa)+mu*beta) ;

                        XK = -NFX + bk*EXPINV(NX,X)  ;
                     
                                X = NX ;
                         
                               FX = NFX ;    
                               
                               
             else
   
                           beta = beta2 ;
                          
%          fprintf('INERTIAL: Norm of (X_1X_3, X_2X_3, X_3^2-1) %d \n',beta)  
                    
                          
                             X = Z ;
                          
                            FX = FZ ;

             end       
                     

                                          
                           k = k + 1 ;
                           
                          

            end              
              
              
                       rdist = DIST(XS, X) ;   
                
                            iterk = k ;    
  
                     TimeCost = etime(clock,t0) ;

                            
          fprintf('\n');
          fprintf('INERTIAL: Number of Iterations %d \n', iterk)
          fprintf('INERTIAL: Number of function estimations %d \n',f_eval)  
          fprintf('INERTIAL: Initial norm of (X_1X_3, X_2X_3, X_3^2-1) ==========%d \n',beta0)
          fprintf('INERTIAL: Final norm of (X_1X_3, X_2X_3, X_3^2-1) ==========%d \n',beta)
          fprintf('INERTIAL: Riemannian distance of X_k and X_* ==========%d \n',rdist)
          fprintf('INERTIAL: Computing time used ========== %d \n', TimeCost)                
                          
         

end

