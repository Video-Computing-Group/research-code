%% Code Author: Ahmed Tashrif Kamal - tashrifahmed@gmail.com
% http://www.ee.ucr.edu/~akamal
% no permission necessary for non-commercial use
% Date: 4/27/2013

%% 
classdef MTIC < handle
    properties
        x;
        x_;
        J;
        J_;
        
        V;
        Vtemp;
        W;
        Wtemp;
        v;
        vtemp;
        
        u;
        U;
        G;
        beta0;
    end
    
    methods
        
        %% constructor
        function a = MTIC(Nc,Nt,p,T,eta,xa,Pinv)
            if (nargin > 0)
                a(Nc,Nt) = MTIC;
                
                
                for i=1:Nc
                    for j=1:Nt
                        a(i,j).x_= xa{j}(:,1) + eta{j};
                        a(i,j).J_ = Pinv;
                        
                        a(i,j).x = zeros(p,T);
                        a(i,j).J = zeros(p,p);
                        a(i,j).v = zeros(p,1);
                        a(i,j).V = zeros(p,p);
                        a(i,j).W = zeros(p,p);
                        a(i,j).vtemp = zeros(p,1);
                        a(i,j).Vtemp = zeros(p,p);
                        a(i,j).Wtemp = zeros(p,p);
                        a(i,j).u = zeros(p,1);
                        a(i,j).U = zeros(p,p);
                        a(i,j).G = zeros(p,p);
                        a(i,j).beta0 = 0;
                    end
                end
            end
        end %end constructor
        
        
        %% data association
        function dataAssociation(a,Nc,Nt,p,m,z,zCount,R,H,FOV,Pg,lamdaf)
            
            for i = 1:Nc
                for j=1:Nt
                    %compute S
                    S = R+H/a(i,j).J_*H';

                    %compute Pd
                    Pd = findPd(H*a(i,j).x_,S,FOV(:,i));
                    
                    %compute b
                    b = 2*pi*lamdaf*sqrt(det(S))*(1-Pd*Pg)/Pd;
                    
                    %initialize stuff
                    y = zeros(m,1);
                    Ptilde = zeros(m,m);
                    
                    if zCount(i) > 0
                        %compute s
                        sum_s = 0;
                        s = zeros(zCount(i),1);
                        for n = 1:zCount(i)
                            ztilde = z{i}(:,n) - H*a(i,j).x_;
                            s(n) = exp(-ztilde' / S * ztilde/2);
                            if s(n)<.0001
                                s(n) = 0;
                            end
                            sum_s = sum_s + s(n);
                        end %for m
                        
                        %compute beta
                        beta = zeros(zCount(i),1);
                        for n=1:zCount(i)
                            beta(n) = s(n)/(b+sum_s);
                            y = y + beta(n)*z{i}(:,n);
                            Ptilde = Ptilde + beta(n)*(z{i}(:,n)-H*a(i,j).x_) * (z{i}(:,n)-H*a(i,j).x_)';
                        end %for m
                        a(i,j).beta0 = b/(b+sum_s);
                        Ptilde = Ptilde - (y - ( 1- a(i,j).beta0 )* H*a(i,j).x_) * (y-( 1- a(i,j).beta0 )*H*a(i,j).x_)';
                        
                        if a(i,j).beta0 < 1
                            % U and u
                            a(i,j).U = H'/R*H;
                            a(i,j).u = H'/R*y;
                            
                            %compute G
                            K = a(i,j).J_\H'/S;
                            C = Ptilde - ( 1- a(i,j).beta0 )*S;
                            a(i,j).G = -a(i,j).J_*K/( inv(C) + K'*a(i,j).J_*K  )*K'*a(i,j).J_;
                        else
                            a(i,j).U = zeros(p,p);
                            a(i,j).u = zeros(p,1);
                            a(i,j).G = zeros(p,p);
                        end
                        
                    else %i.e. if zCount(i) = 0
                        a(i,j).U = zeros(p,p);
                        a(i,j).u = zeros(p,1);
                        a(i,j).G = zeros(p,p);
                    end %end if
                    
                    
                end %end j
            end %end i
            
        end
        
        
        
        
        %% prepData
        function prepData(a,Nc,Nt)
            for i=1:Nc
                for j=1:Nt
                    a(i,j).v = a(i,j).u + (1/Nc * a(i,j).J_ + a(i,j).beta0 *  a(i,j).U) * a(i,j).x_ ;
                    a(i,j).V = 1/Nc * a(i,j).J_ + a(i,j).U;
                    a(i,j).W = 1/Nc * a(i,j).J_ + a(i,j).G;
                end
            end
        end
        
        
        %% consensus
        function consensus(a,Nc,Nt,p,K,eps,E)
            for k=1:K
                for i=1:Nc
                    for j=1:Nt
                        a(i,j).vtemp = zeros(p,1);
                        a(i,j).Vtemp = zeros(p,p);
                        a(i,j).Wtemp = zeros(p,p);                        
                        Delta = sum(E(i,:));
                        for i_=1:Nc
                            if E(i,i_)
                                a(i,j).vtemp = a(i,j).vtemp + a(i_,j).v;
                                a(i,j).Vtemp = a(i,j).Vtemp + a(i_,j).V;
                                a(i,j).Wtemp = a(i,j).Wtemp + a(i_,j).W;
                            end
                        end
                        a(i,j).vtemp =   (1-Delta*eps)*a(i,j).v  + eps*a(i,j).vtemp;
                        a(i,j).Vtemp =   (1-Delta*eps)*a(i,j).V  + eps*a(i,j).Vtemp;
                        a(i,j).Wtemp =   (1-Delta*eps)*a(i,j).W  + eps*a(i,j).Wtemp;
                    end
                end
                
                for i=1:Nc
                    for j= 1:Nt
                        a(i,j).v = a(i,j).vtemp;
                        a(i,j).V = a(i,j).Vtemp;
                        a(i,j).W = a(i,j).Wtemp;
                    end
                end
                
            end
            
        end
        
        %% Estimate
        function estimate(a,Nc,Nt,p,t,Pinv,Phi,Q)
            for i=1:Nc
                for j=1:Nt
                    a(i,j).x(:,t) = a(i,j).V \ a(i,j).v;
                    
                    if all(eig(a(i,j).W)>0) % to make sure the cov is positive
                        a(i,j).J = Nc * a(i,j).W;
                    else
                        a(i,j).J = Pinv;
                    end
                    
                    a(i,j).x_ = Phi * a(i,j).x(:,t);
                    
                    P = Phi / a(i,j).J * Phi' + Q;
                    %fix symmetry
                    P(2,1) = P(1,2);
                    P(3,1) = P(1,3);
                    P(4,1) = P(1,4);
                    P(3,2) = P(2,3);
                    P(4,2) = P(2,4);
                    P(4,3) = P(3,4);                      
                    
                    
                    a(i,j).J_ = eye(p,p)/ ( P );
                    
                  
                end
            end
        end
        
  
        
        
        
        
    end %methods
    
end