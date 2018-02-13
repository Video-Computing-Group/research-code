% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Generates a trunctated version of a monotonic submodular function
% If Ftrunc = sfo_fn_trunc(F,thresh), then
% Ftrunc(A) = min(F(A),thresh)
% Example: See sfo_fn.m and the tutorial script for more information
function F = sfo_fn_trunc(oldF,thresh)
F.oldF = oldF;
F.thresh = thresh;
F = class(F,'sfo_fn_trunc',sfo_fn);
