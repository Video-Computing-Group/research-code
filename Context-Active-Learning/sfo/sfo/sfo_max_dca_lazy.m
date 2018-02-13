% The data-correcting algorithm for maximizing general submodular functions
% of Goldengorin et al.
% Implemented by Andreas Krause (krausea@gmail.com), using lazy evaluations
% This is a worst-case exponential branch & bound algorithm
%
% function A = sfo_max_dca_lazy(F,V,opt)
% F: submodular function
% V: index set
% opt (optional): parameter struct, depending on:
%
% dca_approx: specified accuracy level
% dca_max_card: constraint on #elements; requires F to be monotonic!
%
% returns solution A' s.t. F(A') >= max_A F(A) - alpha
%
% Example: A = sfo_max_dca_lazy(@sfo_fn_example,1:2)

function A = sfo_max_dca_lazy(F,V,opt)
global sfo_max_dca_best_value; 

if ~exist('opt','var')
    opt = sfo_opt;
end
alpha = sfo_opt_get(opt,'dca_approx',0);

sfo_max_dca_best_value = -inf;

if isfield(opt,'dca_max_card')
    k = sfo_opt_get(opt,'dca_max_card');
    [AG,scores] = sfo_greedy_lazy(F,V,k);
    optbound = sfo_maxbound(F,V,AG,k);
    % create new submodular function that's negative if |A|>k
    Fpenalty = sfo_fn_wrapper(@(A) optbound*max(0,length(A)-k));
    FT = sfo_fn_lincomb({F,Fpenalty},[1 -1]);
    sfo_max_dca_best_value = sfo_max_dca_update_best(sfo_max_dca_best_value,scores(end),AG);
else
    k = -1;
    FT = F;
    sfo_max_dca_best_value = F([]);
end

A = sfo_max_dca_rec(FT,[],V,k,alpha,[],[]);


%% The recursive procedure of the Data Correcting Algorithm
function A = sfo_max_dca_rec(F,S,T,k,alpha,boundp,boundm)
global sfo_max_dca_best_value;

[S,T,boundp,boundm] = sfo_max_dca_pp(F,S,T,boundp,boundm);

optbound = sort(boundp,'descend'); optbound = optbound(optbound>=0); 
if k>0 %maximize monotonic function
    optbound = F(S)+sum(optbound(1:min(k,length(optbound))));
else %maximize general function
    optboundm = sort(boundm,'descend'); optboundm = optboundm(optboundm>=0); 
    optboundm = F(T)+sum(optboundm);        
    optboundp = F(S)+sum(optbound);    
    optbound = min(optboundm,optboundp);
end
if optbound<sfo_max_dca_best_value % can break since we already have better solution
    A = S;
    return     
end
if length(S)==length(T)
    A = S;
    return 
end
[kpmax,deltapmax,boundp] = sfo_max_delta_lazy(F,S,T,+1,boundp);
[kmmax,deltammax,boundm] = sfo_max_delta_lazy(F,S,T,-1,boundm);

if deltapmax<=deltammax
    boundp(kpmax) = -inf; boundm(kpmax) = -inf;
    if deltapmax<=alpha
        A = sfo_max_dca_rec(F,S,sfo_setdiff_fast(T,kpmax), k, alpha-deltapmax, boundp, boundm);
        sfo_max_dca_best_value = sfo_max_dca_update_best(sfo_max_dca_best_value,F(A),A);
    else
        A1 = sfo_max_dca_rec(F,[S kpmax],T, k, alpha, boundp, boundm);
        sfo_max_dca_best_value = sfo_max_dca_update_best(sfo_max_dca_best_value,F(A1),A1);
        A2 = sfo_max_dca_rec(F,S,sfo_setdiff_fast(T,kpmax), k, alpha, boundp, boundm);
        sfo_max_dca_best_value = sfo_max_dca_update_best(sfo_max_dca_best_value,F(A2),A2);
        if F(A1)>=F(A2)
            A = A1;
        else
            A = A2;
        end
        return        
    end
else
    boundp(kmmax) = -inf; boundm(kmmax) = -inf;
    if deltammax<=alpha
        A = sfo_max_dca_rec(F,[S kmmax],T, k, alpha-deltammax, boundp, boundm);
        sfo_max_dca_best_value = sfo_max_dca_update_best(sfo_max_dca_best_value,F(A),A);
    else
        A1 = sfo_max_dca_rec(F,[S kmmax],T, k, alpha, boundp, boundm);
        sfo_max_dca_best_value = sfo_max_dca_update_best(sfo_max_dca_best_value,F(A1),A1);
        A2 = sfo_max_dca_rec(F,S,sfo_setdiff_fast(T,kmmax), k, alpha, boundp, boundm);
        sfo_max_dca_best_value = sfo_max_dca_update_best(sfo_max_dca_best_value,F(A2),A2);
        if F(A1)>=F(A2)
            A = A1;
        else
            A = A2;
        end
        return        
    end
    
end

%% Preliminary Preservation algorithm of Goldengorin et al
% extended using lazy evaluations

function [S,T,boundp,boundm] = sfo_max_dca_pp(F,S,T,boundp,boundm)

if length(S)==length(T)
    return
end
while (length(S)~=length(T))
    [kpmax,deltapmax,boundp] = sfo_max_delta_lazy(F,S,T,+1,boundp);
    if deltapmax<=0
        % bounds remain valid, except bounds(kpmax)=-inf
        boundp(kpmax)=-inf;
        boundm(kpmax)=-inf;
        T = sfo_setdiff_fast(T,kpmax); 
    else
        [kmmax,deltammax,boundm] = sfo_max_delta_lazy(F,S,T,-1,boundm);
        if deltammax<=0
            % bounds remains valid, except bounds(kmmax)=-inf
            boundp(kmmax)=-inf;
            boundm(kmmax)=-inf;
            S = [S kmmax]; 
        else
            return
        end        
    end
end

%% update and output stats
function sfo_max_dca_best_value = sfo_max_dca_update_best(sfo_max_dca_best_value,FA,A)
if FA>sfo_max_dca_best_value
    sfo_max_dca_best_value = FA;
    disp(sprintf('new best solution: F(A) = %f. Set = %s',FA,sprintf('%d ',A)));
end
