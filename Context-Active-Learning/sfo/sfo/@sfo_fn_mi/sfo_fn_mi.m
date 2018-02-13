% Computes the Gaussian mutual information between a set and its complement
% Author: Andreas Krause (krausea@gmail.com)
%
% function mi = sfo_fn_mi(sigma,V) 
% sigma: Covariance Matrix
% set: the ground set
%
% Example: F = sfo_fn_mi(0.5*eye(3)+0.5*ones(3),1:3); F(2)

function F = sfo_fn_mi(sigma,V)
F.sigma = sigma;
F.V = V;

F.indsA = [];
F.invAc = [];
F.indsAc = [];
F.cholA = [];

F = class(F,'sfo_fn_mi',sfo_fn);
F = set(F,'current_set',-1);
