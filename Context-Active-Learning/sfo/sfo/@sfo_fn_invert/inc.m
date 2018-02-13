% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
function new_val = inc(F,A,el)
A = sfo_unique_fast(A);
F = init(F,A);
if sum(A==el)>0
    new_val = get(F,'current_val');
    return
end
new_val = dec(F.F, sfo_setdiff_fast(F.V, A), el);
new_val = new_val-F.FV;
