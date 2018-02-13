% The Lovasz extension [Lovasz '83]
% Implementation by Andreas Krause (krausea@gmail.com)
%
% function x = sfo_lovaszext(F,V,w)
% F: Submodular function
% V: index set
% w: weight vector to evaluate Lovasz extension at
%
% Example: x = sfo_lovaszext(@sfo_fn_example,1:2,[0,1])

function x = sfo_lovaszext(F,V,w)
x = w*sfo_polyhedrongreedy(F,V,w)';
