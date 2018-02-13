% Implementation by Andreas Krause (krausea@gmail.com)
% Truncates the marginal variance reduction at variance level c
% Example: See sfo_fn.m and the tutorial script for more information
function F = trunc(F,c)
F = sfo_fn_varred_trunc(F.sigma,F.V,c);

