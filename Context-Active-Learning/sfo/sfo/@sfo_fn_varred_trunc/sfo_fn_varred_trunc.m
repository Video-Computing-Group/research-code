% Implementation by Andreas Krause (krausea@gmail.com)
% Computes the average truncated variance reduction. To be used with
% sfo_saturate.m
% Example: See sfo_fn.m and the tutorial script for more information
function F = sfo_fn_varred_trunc(sigma,V,threshold)
F.sigma = sigma;
F.trunc_thresh = threshold;
F.varPrior = sum(max(diag(sigma)-threshold,0)); %average truncated variance before conditioning
F.V = V;
F.Ainv = [];
F = class(F,'sfo_fn_varred_trunc',sfo_fn);
F = set(F,'current_set',-1);
