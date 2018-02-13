% The polyhedron greedy algorithm [Edmonds '71]
% Implementation by Andreas Krause
%
% function x = sfo_polyhedrongreedy(F,V,w)
% F: Submodular function
% V: index set
% w: weight vector, w(i) is weight of V(i)
%
% Example:
%   x = sfo_polyhedrongreedy(@sfo_fn_example,1:2,sfo_charvector(1:2,1))

function x = sfo_polyhedrongreedy(F,V,w)
n = length(V);
[w,I] = sort(w,'descend');
x = zeros(1,n);
A = [];
Fold = F([]);
for i = 1:n
    Anew = [A V(I(i))];
    x(I(i)) = F(Anew)-Fold;
    A = Anew;
    Fold = Fold + x(I(i));
end
