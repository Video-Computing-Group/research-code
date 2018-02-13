% Computes the Gaussian entropy
% Author: Andreas Krause (krausea@gmail.com)
%
% function H = sfo_fn_entropy(sigma,set) 
% sigma: Covariance Matrix
% set: the subset of rows
%
% Example: F = sfo_fn_entropy(0.5*eye(3)+0.5*ones(3),1:3);

function F = sfo_fn_entropy(sigma,V)
F.sigma = sigma;
F.V = V;

F.indsA = [];
F.cholA = [];

F = class(F,'sfo_fn_entropy',sfo_fn);
F = set(F,'current_set',-1);
