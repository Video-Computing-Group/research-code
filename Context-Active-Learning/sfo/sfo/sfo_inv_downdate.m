% Andreas Krause (krausea@gmail.com)
% Updates a precision matrix by removing a variable
% Based on implementation by Ram Rajagopal
%
% function precNew = sfo_inv_downdate(prec,k)
% prec: old precision (inverse covariance) matrix
% k: index of variable (column) to be removed
%
% Example: See use in @sfo_fn_mi

function precNew = sfo_inv_downdate(prec,k)
n = size(prec,1);
P = [1:(k-1) (k+1):n k];
prec = prec(P,P);

u = -prec(1:(n-1),n);
v = prec(1:(n-1),n)/prec(n,n);
Cinv = prec(1:(n-1),1:(n-1));

precNew = Cinv + u*v';
