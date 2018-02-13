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
        
if (isempty(A))
    sigmaXgA = F.sigma(el,el);
else
	tmpA = (F.cholA'\F.sigma(A,el));
    sigmaXgA = F.sigma(el,el)-F.sigma(el,A)*(F.cholA\tmpA);
end

H = 1/2*log2(2*pi*exp(1)) + 1/2*sum(log2(sigmaXgA));
Hn = log2(sqrt(2*pi*exp(1)*F.noise));

newScore = get(F,'current_val')+H-Hn;
