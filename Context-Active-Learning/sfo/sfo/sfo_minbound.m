% Get bound on suboptimality / certificate of optimality [Edmonds '71]
% Implementation by Andreas Krause (krausea@gmail.com)
%
% function bound = sfo_minbound(F,V,A)
% F: Submodular function
% V: index set
% A: set to be tested
% Returns bound <= min_A F(A)
%
% Example: bound = sfo_minbound(sfo_fn_example,1:2,[1])

function bound = sfo_minbound(F,V,A)
w = sfo_charvector(V,A);
xw = sfo_polyhedrongreedy(F,V,w); %go inside polyhedron
xw = min(xw,0);
bound = sum(xw);
