% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function [F,mi] = init(F,sset)
sset = sfo_unique_fast(sset);
if ~isequal(sset,get(F,'current_set'))
    Ac = sfo_setdiff_fast(F.V,sset);
    F.invAc = inv(F.sigma(Ac,Ac));
    F.cholA = chol(F.sigma(sset,sset));
    F.indsA = sset;
    F.indsAc = Ac;
    
    if isempty(sset)
        mi = 0;
    else
        sigmaA = F.sigma(sset,sset)+(1e-10)*(eye(length(sset)));
        sigmaAcomp = F.sigma(Ac,Ac)+(1e-10)*(eye(length(Ac)));
        sigmaAAcomp = F.sigma(sset,Ac);

        sigmaAcond = sigmaA-sigmaAAcomp * (sigmaAcomp \ sigmaAAcomp');
        H = 1/2*log2((2*pi*exp(1))^size(sigmaA,1)) + 1/2*sfo_logdet(sigmaA);
        Hcond = 1/2*log2((2*pi*exp(1))^size(sigmaAcond,1)) + 1/2*sfo_logdet(sigmaAcond);
        mi = H-Hcond;
        
    end
    F = set(F,'current_val',mi,'current_set',sset);
else
    mi = get(F,'current_val');
end
