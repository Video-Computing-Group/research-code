% Implementation by Andreas Krause (krausea@gmail.com)
% Variance reduction in Gaussian linear models
% sigma is the covariance matrix
% V is the ground set
% computes the expected mean squared prediction error (trace of posterior
% covariance)
% Supports the method trunc to be used in conjunction with sfo_saturate
% Example: See sfo_fn.m and the tutorial script for more information
function F = sfo_fn_varred(sigma,V)

if nargin==0
    F.sigma = 0;
    F.V = 0;
    F.Ainv = [];
    F.AAc = [];
    F.trD = [];
else
    F.sigma = sigma;
    F.V = V;
    F.Ainv = [];
    F.AAc = [];
    F.trD = [];
end

F = class(F,'sfo_fn_varred',sfo_fn);
F = set(F,'current_set',-1);

