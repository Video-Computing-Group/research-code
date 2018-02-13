% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function [F,v] = init(F,A)
A = sfo_unique_fast(A);
v = 0;
for i = 1:length(F.Fs)
    F.Fs{i} = init(F.Fs{i},A);
    v = v+F.weights(i)*get(F.Fs{i},'current_val');
end
F = set(F,'current_set',A,'current_val',v);
