% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
function [F,val] = init(F,sset)
sset = sfo_unique_fast(sset);
if ~isequal(sset,get(F,'current_set'))

    if (isempty(sset))
        val = 0;
    else

        comp = sfo_setdiff_fast(F.V,sset);
        F.Ainv = inv(F.sigma(sset,sset));

        varPost = 0;

        for i = 1:length(comp)
            s = comp(i);
            v = F.sigma(s,s)-F.sigma(s,sset)*F.Ainv*F.sigma(sset,s);
            varPost = varPost + max(v-F.trunc_thresh,0);
        end       
        val = (F.varPrior-varPost)/length(F.V);
    end
    F = set(F,'current_set',sset,'current_val',val);
else
    val = get(F,'current_val');
end


