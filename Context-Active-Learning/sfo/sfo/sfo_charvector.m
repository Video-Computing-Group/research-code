% returns characteristic vector of A
% Author: Andreas Krause (krausea@gmail.com)
%
% function w = sfo_charvector(V,A)
% V: index set
% A: subset of V
%
% Example: sfo_charvector([2,1,4],[4,2])=[1 0 1]

function w = sfo_charvector(V,A)
w = ismember(V,A);
