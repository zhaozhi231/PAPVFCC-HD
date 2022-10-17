function  [ PW ] = PC(W)


%                Y = randn(2,1) ;     
%                 
%             Y = [ Y ; sqrt(Y'*Y+1) ] ;
%             
%             LB = [-Inf ; -Inf ; 1];
%             
%             UB = [ Inf ;  Inf ; 2];  
%             
% 
%          options = optimoptions('fmincon','Algorithm','sqp'); % run interior-point algorithm
%          PX = fmincon(@(X) myfun(X,W),Y,[],[],[],[],LB,UB,@(X) mycon(X),options)  ;
%             
%   
%          function f = myfun(X,W)       
%          
%              f = acosh(max(-LORENTZ(X,W),1)) ;          
%          
%          end
%      
%      
%          function [c,ceq] = mycon(X)
%         
%               c =  [] ;
%                             
%            ceq = LORENTZ(X,X)+1 ; 
%            
%   
%          end
     

%                W = randn(2,1) ;     
%                 
%             W = [ W ; sqrt(W'*W+1) ] ;
%             
%              aa =  LORENTZ(W,W) ;

    if LORENTZ(W,W) ~=-1 || W(3)<1 || W(3)>2

               Y = randn(2,1) ;     
                
            Y = [ Y ; sqrt(Y'*Y+1) ] ;
            
              bb =  LORENTZ(W,W) ;

            
            LB = [-Inf ; -Inf ; 1];
            
            UB = [ Inf ;  Inf ; 2];  
            

         options = optimoptions('fmincon','Display','off','Algorithm','sqp'); % run interior-point algorithm
         PW = fmincon(@(X) myfun(X,W),W,[],[],[],[],LB,UB,@(X) mycon(X),options)  ;
            
%           cc =  LORENTZ(PX,PX) ;


    else
            
        
            PW = W ;
            
    end
         
         function f = myfun(X,W)      
 
             
%              f = acosh(max(-LORENTZ(X,W),1))  ;      
         
             f = -LORENTZ(X,W) ;        
         
         end
     
     
         function [c,ceq] = mycon(X)
        
              c =  [] ;
                            
           ceq = LORENTZ(X,X)+1 ;
           
  
         end
         
         




%                 Y = randn(2,1) ;   
%                 
%          options = optimoptions('fmincon','Algorithm','sqp'); % run interior-point algorithm
%          SX = fmincon(@(X) myfun(X,W),Y,[],[],[],[],[],[],@(X) mycon(X),options)  ;
%          
%              PX = [ SX ; sqrt(SX'*SX+1) ] ;
%             
%   
%          function f = myfun(X,W)       
%          
%              f = sqrt(1+X(1)^2+X(2)^2)*W(3) - X(1)*W(1) - X(2)*W(2) ;          
%          
%          end
%      
%      
%          function [c,ceq] = mycon(X)
%     
%               c =  [ 4 - (1+X(1)^2+X(2)^2) (1+X(1)^2+X(2)^2) -9 ] ;
%                             
%              ceq = [ ] ; 
%            
%   
%          end
         
 end
