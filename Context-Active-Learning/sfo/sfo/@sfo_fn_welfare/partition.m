% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
% groups elements of the same color in a bucket
function A_part = partition(F,As)
m = length(F.Fs);
groups = mod(As,m)+1;
els = floor(As/m);
A_part = {};
for i = 1:m
    A_part{i} = els(groups==i);
end
