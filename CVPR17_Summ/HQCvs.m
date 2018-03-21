% -------------------------------------------------------------------------
% Main Program for Collaborative video summarization
% Input: set of feature matrices: X and X_n
% Similarity matrices: D and D_n
% Regularization parameters: alpha, lambda_s, lambda_d, beta, maxIter
% Output: Sparse coefficient matrix: Z_c
% Copyright @ Rameswar Panda, 2016
% Ref: Efficient and Robust Feature Selection via Joint L2,1-Norms
% Minimization, NIPS'10
% -------------------------------------------------------------------------
function Zc = HQCvs(X,X_n,D,D_n,alpha,lambda_s,lambda_d,beta,maxIter)


n = size(X,2);
p = ones(n,1); % initialization
q = ones(n,1);
r = ones(n,1);

for i = 1:maxIter
    
    P = diag(p);
    Q = diag(q);
    R = diag(r);
    
    Z1 = (X'*X + 2*lambda_s*P + 2*beta*R)\(X'*X - lambda_d*D);
    Z2 = (alpha*X'*X + 2*lambda_s*Q + 2*beta*R)\(alpha*X'*X_n - lambda_d*D_n);
    Zc = horzcat(Z1,Z2);
    
    pt = sqrt(sum(Z1.*Z1,2)+0.0001);
    p = 1./(2.*pt);
   
    qt = sqrt(sum(Z2.*Z2,2)+0.0001);
    q = 1./(2.*qt);
    
    rt = sqrt(sum(Zc.*Zc,2)+0.0001);
    r = 1./(2.*rt);       
end
end



    
    

