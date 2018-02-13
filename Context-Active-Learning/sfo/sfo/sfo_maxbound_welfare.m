% Getting an online bound on the optimal solution for submodular welfare
% Implementation by Andreas Krause (krausea@gmail.com)
%
% function bound = sfo_maxbound_welfare(F,V,As)
% F: Submodular function
% V: index set
% As: current sets (cell array of sets)
% Returns: bound on optimal solution
%
% Example: 
%   F = sfo_fn_entropy(ones(5)+eye(5),1:5); 
%   bound = sfo_maxbound_welfare(F,1:5,{[1,2],[3 4 5]})

function bound = sfo_maxbound_welfare(F,V,As)
m = length(As);
k = length([As{:}]);
n = length(V);

delta = zeros(n,m);
curVal = 0;

for i = 1:m
    F = init(F,As{i});
    valI = F(As{i});
    
    for s = 1:n
        delta(s,i)=inc(F,As{i},V(s))-valI;
    end
    curVal = curVal+valI;
end
if m>1
    delta = max(delta');
end
deltas = sort(delta,'descend');
bound = curVal/m + sum(deltas(1:k))/m;
