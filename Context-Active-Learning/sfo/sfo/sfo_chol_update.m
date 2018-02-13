% Andreas Krause (krausea@gmail.com)
% This function is used to update the Cholesky decomposition of the
% X'X(A,A) matrix when a new variable is added to the model
% It checks whether the addition causes X'X(A,A) to become singular and
% returns a flag indicating that
% Based on implementation by Ram Rajagopal
%
% function [newChol,rankok] = sfo_chol_update(sigma, A, new_indices, oldChol, my_eps) 
% sigma: joint covariance matrix
% A: indices into sigma of old variables (representing oldChol)
% new_indices: indices into sigma that are added to oldChol
% oldChol: cholesky decomposition if sigma(A,A)
% my_eps: tolerance for checking rank
%
% Example: See use in @sfo_fn_entropy

function [newChol,rankok] = sfo_chol_update(sigma, A, new_indices, oldChol, my_eps)
xtR = sigma(A,new_indices);
xtx = sigma(new_indices,new_indices);

if(isempty(oldChol))
 rankok = 1;
 newChol   = chol(xtx);
 return;
end;

r   = (oldChol')\xtR;
rpp = xtx - r'*r;
if(rpp <= my_eps)
 rankok = 0;
 newChol = [];
else
 rankok = 1;
 newChol   = [oldChol r;zeros(1, size(oldChol,2)) sqrt(rpp)];
end
