% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function newScore = inc(F,A,el)
A = sfo_unique_fast(A);
F = init(F,A);

if sum(A==el)>0
    newScore = get(F,'current_val');
    return;
end

Ac = sfo_setdiff_fast(F.V,[A el]);
        
pos = find(F.indsAc == el);
invAc = sfo_inv_downdate(F.invAc, pos);

if (isempty(A))
    sigmaXgA = F.sigma(el,el);
else
    sigmaXgA = F.sigma(el,el)-F.sigma(el,A)*(F.cholA\(F.cholA'\F.sigma(A,el)));
end

sigmaXgAc = F.sigma(el,el)-F.sigma(el,Ac)*(invAc*F.sigma(Ac,el));

H = 1/2*log2(2*pi*exp(1)) + 1/2*sum(log2(sigmaXgA));
Hcond = 1/2*log2(2*pi*exp(1)) + 1/2*sum(log2(sigmaXgAc));

mi = H-Hcond;

newScore = get(F,'current_val')+mi;
