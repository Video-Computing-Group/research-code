% implements Queyranne's algorithm for minimizing symmetric submodular
% functions
% author: Andreas Krause (krausea@gmail.com)
%
% function R = sfo_queyranne(F,V)
% F is the symmetric submodular function
% V is the index set
% Returns an optimal solution to min F(A) s.t. 0<|A|<n
%
% Example: 
%   G = [1 1 0; 1 0 1; 0 1 1]; F = sfo_fn_cutfun(G);
%   A = sfo_queyranne(F,1:3)

function R = sfo_queyranne(F,V)
n = length(V);
S = {};
for j = 1:n
    S{j}=V(j);
end
s = zeros(1,n-1);
A = {};
inew = 1:n; 
for h = 1:(n-1)
    Fnew = @(A) F([S{A}]);
    % find a pendant pair
    [t,u]=sfo_pendentpair(Fnew,inew);
    
    % this gives a candidate solution
    A{h} = S{u};
    s(h) = Fnew(u);
    S{t}=[S{t} S{u}];
    inew = sfo_setdiff_fast(inew, u);
    S{u}= -S{u};
end
% return best candidate solution
i = find(s==min(s),1);
R = A{i};

%% 
% implements the pendant pair finding subroutine of queyranne's algorithm
% (Queyranne '95)
% F is the submodular function
% inds is an array of indices; (typically, 1:n)

function [s,t] = sfo_pendentpair(F,V)
x = 1; % x is a starting element, pointing into V;
vnew = V(x);
n = length(V);
Wi = [];
used = zeros(1,n);
used(x) = 1;
for i = 1:(n-1)
    vold = vnew;
    Wi = [Wi vold];
    % now update the keys
    keys = 1e99*ones(1,n);
    for j = 1:n
        if used(j)
            continue;
        end
        keys(j) = F([Wi V(j)]) - F(V(j));
    end
    argmin = find(keys==min(keys),1);
    vnew = V(argmin);
    used(argmin)=1;
end
s = vold;
t = vnew;
