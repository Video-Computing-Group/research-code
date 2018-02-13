% computes the Gaussian maximum posterior variance after conditioning on
% subset set
% Author: Andreas Krause (krausea@gmail.com)
%
% function v = sfo_eval_maxvar(sigma,set)
% sigma: Covariance matrix
% set: subset of variables
%
% Example: v = sfo_eval_maxvar(0.5*eye(3)+0.5*ones(3),1);

function v = sfo_eval_maxvar(sigma,set)
n=size(sigma,1);
set = sfo_unique_fast(set);

if (isempty(set))
    v = max(diag(sigma));
else    
    comp = setdiff(1:n,set);

    sigmaA = sigma(set,set);
    sigmaAcomp = sigma(comp,comp);
    sigmaAcompA = sigma(comp,set);

    sigmaAcompcond = sigmaAcomp-sigmaAcompA * (sigmaA \ sigmaAcompA');
    v = max(diag(sigmaAcompcond));
end
