% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
%
% function F = sfo_fn_detect(detmat,V)
%
% detmat is a N x D (sparse) matrix where N is #sensors, D is #scenarios;
% detmat(i,j) is benefit if sensor i detects scenario j
function F = sfo_fn_detect(detmat,V)

if nargin==0
    F.detmat = 0;
    F.V = 0;
    F.curmax=0;
    F.marginals=0;
else
    F.detmat = detmat';
    F.V = V;
    F.curmax = zeros(size(F.detmat,1),1);
    F.marginals = sum(F.detmat,1);
end

F = class(F,'sfo_fn_detect',sfo_fn);
F = set(F,'current_set',-1);
