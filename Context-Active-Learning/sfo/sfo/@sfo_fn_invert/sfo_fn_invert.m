% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Given a submodular function G, this function represents the "inverse"
% F(A) = G(V\A)
% Example: See sfo_fn.m and the tutorial script for more information
function F = sfo_fn_invert(oldF,V)
F.F = init(oldF,V);
F.V = V;
F.FV = get(F.F,'current_val');
F = class(F,'sfo_fn_invert',sfo_fn);
