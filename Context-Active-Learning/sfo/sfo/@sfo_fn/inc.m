% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function new_val = inc(F,A,el)
F = init(F,[A, el]);
new_val = get(F,'current_val');
