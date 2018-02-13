% Andreas Krause (krausea@gmail.com)
% Deletes a variable from the X'X matrix in a Cholesky factorisation R'R =
% X'X. Returns the downdated R. This function is just a stripped version of
% Matlab's qrdelete. 
% Based on implementation by Ram Rajagopal, originally from Kevin Murphy
%
% function R = sfo_chol_downdate(R,j)
% R: old cholesky decomposition
% j: column to be removed
%
% Example: See use in @sfo_fn_mi

function R = sfo_chol_downdate(R,j)
R(:,j) = []; % remove column j
n = size(R,2);
for k = j:n
  p = k:k+1;
  [G,R(p,k)] = planerot(R(p,k)); % remove extra element in column
  if k < n
    R(p,k+1:n) = G*R(p,k+1:n); % adjust rest of row
  end
end
R(end,:) = []; % remove zero'ed out row

