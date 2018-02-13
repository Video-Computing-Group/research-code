% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Creates a (positive) linear combination of submodular functions
% Example: See sfo_fn.m and the tutorial script for more information

function F = sfo_fn_lincomb(Fs,weights)
F.Fs = Fs;
F.weights = weights;
F = class(F,'sfo_fn_lincomb',sfo_fn);
