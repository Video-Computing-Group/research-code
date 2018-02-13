% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
function [F,val] = init(F,sset)
sset = sfo_unique_fast(sset);
if ~isequal(sset,get(F,'current_set'))
    Ac = sfo_setdiff_fast(F.V,sset);
    F.Ainv = inv(F.sigma(sset,sset));
    F.AAc = F.sigma(sset,Ac)*F.sigma(Ac,sset);
    F.trD = trace(F.AAc*F.Ainv);
    val = (trace(F.sigma(sset,sset))+F.trD)/length(F.V); %sfo_fn_sat_varred(F.sigma,sset);
    F = set(F,'current_set',sset,'current_val',val);
else
    val = get(F,'current_val');
end