% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function new_val = inc(F,A,el)
A = sfo_unique_fast(A);
if sum(A==el)>0
    new_val = get(F,'current_val');
    return
end

new_val = 0;
for i = 1:length(F.Fs)
    v = inc(F.Fs{i},A ,el);
    new_val = new_val+F.weights(i)*v;
end
