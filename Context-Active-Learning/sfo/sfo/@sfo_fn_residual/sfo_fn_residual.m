% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Creates a residual submodular function, with the property that if
% Fresid = sfo_fn_residual(F,A), then Fresid(B) = F([A B])-F(A)
% Example: See sfo_fn.m and the tutorial script for more information
function F = sfo_fn_residual(oldF,sset)
sset = sfo_unique_fast(sset);

F.oldF = oldF;
F.sset = sset;
F.ssetVal = oldF(sset);
F = class(F,'sfo_fn_residual',sfo_fn);
