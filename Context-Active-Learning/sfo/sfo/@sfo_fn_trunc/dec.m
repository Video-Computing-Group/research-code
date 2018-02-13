% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
function new_val = dec(F,A,el)
A = sfo_unique_fast(A);
if sum(A==el)==0
    new_val = get(F,'current_val');
    return
end

new_val = dec(F.oldF,A,el);

if F.thresh>=0
    new_val = min(F.thresh,new_val);
end
