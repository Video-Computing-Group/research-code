% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
function scoreNew = inc(F,A,el)
A = sfo_unique_fast(A);
F = init(F,A);

n=length(F.V);
oldScore = get(F,'current_val');

if sum(A == el)>0
    scoreNew = oldScore;
    return
end

Ac = sfo_setdiff_fast(F.V,[A el]);

b = F.sigma(A,el);
newAAc11 = F.AAc - b*b';
newAAc12 = F.sigma(A,Ac)*F.sigma(Ac,el);
newAAc22 = F.sigma(el,Ac)*F.sigma(Ac,el);
newAAc = [newAAc11 newAAc12; newAAc12' newAAc22];

newAinv = sfo_inv_update( F.Ainv, F.sigma(A,el), F.sigma(el,el));

trDNew = 0;
for i = 1:(length(A)+1)
    trDNew = trDNew+newAAc(i,:)*newAinv(:,i);
end

imp = (F.sigma(el,el)+trDNew - F.trD)/n;
scoreNew = oldScore + imp;
