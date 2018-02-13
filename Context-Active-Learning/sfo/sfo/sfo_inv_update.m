% Andreas Krause (krausea@gmail.com)
% Updates a precision matrix by adding a new variable
% Based on implementation by Ram Rajagopal
%
% function precNew = sfo_inv_update(prec,sigma_Ax,sigma_xx)
% prec: old precision (inverse covariance) matrix
% sigma_Ax: column vector of cross-covariances
% sigma_xx: variance of new element
%
% Example: See use in @sfo_fn_mi

function precNew = sfo_inv_update(prec,sigma_Ax,sigma_xx)
r = prec*sigma_Ax;
C2 = sigma_xx - sigma_Ax'*r;

P11 = prec + r*r'/C2;
P12 = -r/C2;
P22 = 1/C2;

precNew = [P11 P12; P12' P22];
