% Implementation by Andreas Krause (krausea@gmail.com)
% To be used by sfo_greedy_welfare
% Example: See sfo_fn.m and the tutorial script for more information
function F = sfo_fn_welfare(Fs)
F.Fs = Fs;
F = class(F,'sfo_fn_welfare',sfo_fn);

