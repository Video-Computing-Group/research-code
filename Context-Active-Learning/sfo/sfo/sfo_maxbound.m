% Getting an online bound on the optimal solution for budgeted maximization
% Implementation by Andreas Krause (krausea@gmail.com)
%
% function bound = sfo_maxbound(F,V,A,B,C)
% F: Submodular function
% V: index set
% A: current set (optional)
% B: Budget
% C: Cost (optional). C(i) is cost of V(i)
% Returns: bound on optimal solution
%
% Example: 
%   F = sfo_fn_entropy(ones(5)+eye(5),1:5); 
%   bound = sfo_maxbound(F,1:5,[1,2],2)

function bound = sfo_maxbound(F,V,A,B,C)
n = length(V);
if isstruct(F)
    F = F.F;
end
if ~exist('A','var')
    % return unconstrained bound
    F0 = F([]);
    boundAdd = F0;
    for i = 1:length(V)
        boundAdd = boundAdd+(F(V(i))-F0);
    end
    FV = F(V);
    boundRem = FV;
    for i = 1:length(V)
        boundRem = boundRem+max(0,(F(V([1:(i-1) (i+1):n]))-FV));
    end
    bound = min(boundAdd,boundRem);
else
    if ~exist('C','var')
        C = ones(1,n);
    end
    deltas = -inf*ones(1,n);
    curVal = F(A);
    for i = 1:n
        if sum(V(i)==A)>0
            continue
        end
        deltas(i) = (F([A i])-curVal)/C(i);
    end
    [deltas,I] = sort(deltas,'descend');
    C = C(I);
    naff = find(cumsum(C)<=B,1,'last'); 
    bound = sum(deltas(1:naff).*C(1:naff));

    frac = (B-sum(C(1:naff)));
    bound = bound+deltas(naff+1)*frac + curVal;
end