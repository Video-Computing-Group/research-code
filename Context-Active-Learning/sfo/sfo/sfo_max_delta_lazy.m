% Helper routine for computing lazy increments / decrements
% Implemented by Andreas Krause (krausea@gmail.com)
%
% function A = sfo_max_delta_lazy(F,S,T,pm,bound)
% F: submodular function
% S,T: sets, S subset of T
% pm: if +1, will compute increments for adding, -1 for removing elements
% bound: upper bound for lazy evaluations. Can be set to [] for
% initialization
%
% Returns best element to add / remove, the increment value + new bounds
%
% Example: argmax = sfo_max_delta_lazy(sfo_fn_example,[],1:2,1,[])

function [argmax,maxval,bound] = sfo_max_delta_lazy(F,S,T,pm,bound)

D = sfo_setdiff_fast(T,S); 
M = max(D(:));
if(isempty(bound))
    bound = -inf*ones(1,M);
    bound(D)=inf;
end

bestimprov = -inf;
if pm == 1
    F = init(F,S);
    FS = get(F,'current_val');
elseif pm == -1
    F = init(F,T);
    FT = get(F,'current_val');
end    

[tmp,order] = sort(bound,'descend');
for test = order
    if bound(test)==-inf
        break;
    end
    if bound(test)>bestimprov
        if pm == 1
            improv = inc(F,S,test)-FS; 
        elseif pm == -1
            improv = dec(F,T,test)-FT;
        end
        bound(test)= improv;
        bestimprov = improv;
    end
end
maxval = max(bound);
argmax = find(bound==maxval,1); % find best delta    
