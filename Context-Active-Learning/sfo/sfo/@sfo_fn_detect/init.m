% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
function [F,v] = init(F,sset)
sset = sfo_unique_fast(sset);
if length(sset)==1
    F = set(F,'current_val',F.marginals(sset)); 
    F.curmax = F.detmat(:,sset);        
elseif ~isequal(sset,get(F,'current_set'))
    if isempty(sset)
        F.curmax = zeros(size(F.detmat,1),1);
    else
        F.curmax = max(F.detmat(:,sset),[],2);        
    end
    F = set(F,'current_val',full(sum(F.curmax)),'current_set',sset);
end
v = get(F,'current_val');
