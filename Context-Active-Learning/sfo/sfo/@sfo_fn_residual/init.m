% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
function [F,v] = init(F,A)
A = sfo_unique_fast(A);
F.oldF = init(F.oldF,sfo_unique_fast([A F.sset]));
v = get(F.oldF,'current_val')-F.ssetVal;
F = set(F,'current_set',A,'current_val',v);
