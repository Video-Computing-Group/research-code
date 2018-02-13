% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
function [F,v] = init(F,A)
A = sfo_unique_fast(A);
F.oldF = init(F.oldF,A);
v = get(F.oldF,'current_val');
if F.thresh>=0
    v = min(F.thresh,v);
end
F = set(F,'current_set',A,'current_val',v);
