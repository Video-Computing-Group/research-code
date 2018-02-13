% Helper for quickly computing the set difference
% Author: Andreas Krause (krausea@gmail.com)
%
% function C=sfo_setdiff_fast(A,B)
% A, B: sets (arrays) of positive integers
%
% Example: C = sfo_setdiff_fast([1 3 7 8], [3 7])

function C=sfo_setdiff_fast(A,B)
mx = max([max(A),max(B)]);
vals = zeros(1,mx);
vals(A)=1;
vals(B)=0;
C = find(vals);
