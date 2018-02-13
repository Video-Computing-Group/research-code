% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
function [F,v] = init(F,A)
A = sfo_unique_fast(A);
F.F = init(F.F,sfo_setdiff_fast(F.V,A));
v =  get(F.F,'current_val')-F.FV;
F = set(F,'current_set',A,'current_val',v);
