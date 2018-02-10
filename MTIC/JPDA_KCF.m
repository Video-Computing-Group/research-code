%% Code Author: Ahmed Tashrif Kamal - tashrifahmed@gmail.com
% http://www.ee.ucr.edu/~akamal
% no permission necessary for non-commercial use
% Date: 4/27/2013

%% 
classdef JPDA_KCF < handle
    properties
        x;
        x_;
        P_;
        M
        
        u;
        U;
        beta0;
        K
        Ptilde;
        b;
        B;
        gamma2;
        sum_dx;
    end
    
    methods
        
        %% constructor
        function a = JPDA_KCF(Nc,Nt,p,m,T,eta,xa,P)
            if (nargin > 0)
                a(Nc,Nt) = JPDA_KCF;
                
                
                for i=1:Nc
                    for j=1:Nt
                        a(i,j).x_= xa{j}(:,1) + eta{j};
                        a(i,j).P_ = P;
                        
                        a(i,j).x = zeros(p,T);
                        a(i,j).M = zeros(p,p);
                        
                        a(i,j).u = zeros(p,1);
                        a(i,j).U = zeros(p,p);
                        a(i,j).beta0 = 0;
                        a(i,j).K = zeros(p,m);
                        a(i,j).Ptilde = zeros(m,m);
                        a(i,j).b = zeros(p,1);
                        a(i,j).B = zeros(p,p);
                        a(i,j).gamma2 = zeros(p,p);
                        a(i,j).sum_dx = zeros(p,1);
                    end
                end
            end
        end %end constructor
        
        
        %% data association
        function dataAssociation(a,Nc,Nt,p,m,z,zCount,R,H,FOV,Pg,lamdaf)
            
            for i = 1:Nc
                for j=1:Nt
                    
                    %assuming in the beginning that none are associated
                    a(i,j).beta0 = 1;
                    
                    %compute S
                    S = R+H*a(i,j).P_*H';
                    
                    %compute Pd
                    Pd = findPd(H*a(i,j).x_,S,FOV(:,i));
                    
                    %compute b
                    b_ = 2*pi*lamdaf*sqrt(det(S))*(1-Pd*Pg)/Pd;
                    
                    %initialize stuff
                    y = zeros(m,1);
                    a(i,j).Ptilde = zeros(m,m);
                    
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
                            beta(n) = s(n)/(b_+sum_s);
                            y = y + beta(n)*z{i}(:,n);
                            a(i,j).Ptilde = a(i,j).Ptilde + beta(n)*(z{i}(:,n)-H*a(i,j).x_) * (z{i}(:,n)-H*a(i,j).x_)';
                        end %for m
                        a(i,j).beta0 = b_/(b_+sum_s);
                        a(i,j).Ptilde = a(i,j).Ptilde - (y - ( 1- a(i,j).beta0 )* H*a(i,j).x_) * (y-(1-a(i,j).beta0 )*H*a(i,j).x_)'; % edited the actual jpda-kcf
                        %a(i,j).Ptilde = a(i,j).Ptilde - (y -  H*a(i,j).x_) * (y-H*a(i,j).x_)'; %just testing the actual jpda-kcf
                        
                        if a(i,j).beta0 < 1
                            % U and u
                            a(i,j).U = H'/R*H;
                            a(i,j).u = H'/R*y;
                            
                            %compute G
                            a(i,j).K = a(i,j).P_*H'/S;
                        else
                            a(i,j).U = zeros(p,p);
                            a(i,j).u = zeros(p,1);
                            a(i,j).K = zeros(p,m);
                        end
                        
                    else %i.e. if zCount(i) = 0
                        a(i,j).U = zeros(p,p);
                        a(i,j).u = zeros(p,1);
                        a(i,j).K = zeros(p,m);
                    end %end if
                    
                    
                end %end j
            end %end i
            
        end
        
        
        
        
        %% commOnce
        function commOnce(a,Nc,Nt,p,E)
            for i=1:Nc
                for j=1:Nt
                    a(i,j).B = a(i,j).U;
                    a(i,j).b = a(i,j).u;
                    for i_ = 1:Nc %if this for loop is not executed, we found that performance increases in JPDA-KCF
                        if E(i,i_) == 1
                            a(i,j).B = a(i,j).B + a(i_,j).U;
                            a(i,j).b = a(i,j).b + a(i_,j).u;
                        end
                    end
                    if a(i,j).beta0<1
                        a(i,j).M = a(i,j).beta0*a(i,j).P_ + (1-a(i,j).beta0) * eye(p,p)/ ( inv(a(i,j).P_) + a(i,j).B ) + a(i,j).K*a(i,j).Ptilde*a(i,j).K';
                    else
                        a(i,j).M = a(i,j).P_;
                    end
                end
            end
        end
        
        
        
        %% Estimate
        function estimate(a,Nc,Nt,p,t,eps,E)
            for i=1:Nc
                for j=1:Nt
                    
                    a(i,j).sum_dx = zeros(p,1);
                    a(i,j).gamma2 = eps * a(i,j).M / ( 1 + sqrt(trace(a(i,j).M'*a(i,j).M)) );
                    
                    for i_ = 1:Nc
                        if E(i,i_) == 1
                            a(i,j).sum_dx = a(i,j).sum_dx + a(i_,j).x_ - a(i,j).x_;
                        end
                    end
                    a(i,j).x(:,t) = a(i,j).x_ +  (inv(a(i,j).P_)+a(i,j).B) \ ( a(i,j).b - (1-a(i,j).beta0)*a(i,j).B * a(i,j).x_ ) + a(i,j).gamma2 * a(i,j).sum_dx;
                end
                
            end
        end

        
        
        %% consensus
        function consensus(a,Nc,Nt,p,K,E,t)
            for k=2:K
                for i = 1:Nc
                    for j=1:Nt
                        a(i,j).sum_dx = zeros(p,1);
                        for i_ = 1:Nc
                            if E(i,i_) == 1
                                a(i,j).sum_dx = a(i,j).sum_dx + a(i_,j).x(:,t) - a(i,j).x(:,t);
                            end
                        end
                    end
                end
                
                for i=1:Nc
                    for j=1:Nt
                        a(i,j).x(:,t)  = a(i,j).x(:,t)  + a(i,j).gamma2 * a(i,j).sum_dx;
                    end
                end
                
            end
            
        end
        
        %% Predict
        function predict(a,Nc,Nt,t,P,Phi,Q)
            for i=1:Nc
                for j=1:Nt
                    % check for positive definiteness
                    if all(eig(a(i,j).M)>0) % to make sure the cov is positive
                        a(i,j).P_ = Phi*a(i,j).M*Phi'+Q;
                    else
                        a(i,j).P_ = P;
                    end
                    
                    %preserve symmetry
                    a(i,j).P_(2,1) = a(i,j).P_(1,2);
                    a(i,j).P_(3,1) = a(i,j).P_(1,3);
                    a(i,j).P_(4,1) = a(i,j).P_(1,4);
                    a(i,j).P_(3,2) = a(i,j).P_(2,3);
                    a(i,j).P_(4,2) = a(i,j).P_(2,4);
                    a(i,j).P_(4,3) = a(i,j).P_(3,4);
                    
                    a(i,j).x_ = Phi*a(i,j).x(:,t);
                    
                    
                end
            end
        end
        
        
        
        
        
        
        
    end %methods
    
end