% Compute the log determinant of a positive definite matrix using stable
% cholesky decomposition
% Author: Andreas Krause (krausea@gmail.com)
%
% function D = sfo_logdet(sigma)
% sigma: Covariance matrix
% set: subset of variables
%
% Example: D = sfo_logdet(0.5*ones(3)+0.5*eye(3))

function D = sfo_logdet(sigma)
L = chol(sigma);
D = 2*sum(log2(diag(L)));
