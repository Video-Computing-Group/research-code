% Computes the Gaussian entropy
% Author: Andreas Krause (krausea@gmail.com)
%
% function H = sfo_fn_entropy(sigma,set) 
% sigma: Covariance Matrix
% set: the subset of rows
%
% Example: F = sfo_fn_mi(0.5*eye(3)+0.5*ones(3));

function F = sfo_fn_infogain(sigma,V,noise)
F.sigma = sigma+eye(length(V))*noise;
F.V = V;
F.noise = noise;

F.indsA = [];
F.cholA = [];

F = class(F,'sfo_fn_infogain',sfo_fn);
F = set(F,'current_set',-1);
