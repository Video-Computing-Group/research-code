% Andreas Krause (krausea@gmail.com)
% implements the (deterministic) local search procedure for maximizing
% nonnegative submodular functions by Feige, Mirrokni, Vondrak
% Obtains a 1/3 approximation to max_A F(A)
% The algorithm implements lazy evaluations for efficiency
%
% function sset = sfo_ls_lazy(F,V,B,C,useCB) 
% F: submodular function
% V: index set
% opt: optional parameter struct, with following entries:
%
% ls_tolerance: required improvement as stopping criterion (default = 0)
% ls_initial_set: starting set;
%
% returns solution sset 
%
% Example: A = sfo_ls_lazy(F_cut_dir,V_G)

function sset = sfo_ls_lazy(F,V,opt) 
if ~exist('opt','var')
    opt = sfo_opt;
end

n=length(V);
TOL = sfo_opt_get(opt,'ls_tolerance',1e-6);
initA = sfo_opt_get(opt,'ls_initial_set',[]);

while true
    % upward pass of greedily adding elements
    F = init(F,[]);
    initA = sfo_greedy_lazy(F,V,inf,sfo_opt({'greedy_initial_sset',initA}));
    cur_val = F(initA);
    if sfo_opt_get(opt,'verbosity_level',0)>0
        fprintf('Value after upward pass: \t%f\n',cur_val);
    end
    
    % downward pass of greedily removing elements
    Finv = sfo_fn_invert(F,initA);
    decA = sfo_greedy_lazy(Finv,initA,inf);
    initA = sfo_setdiff_fast(initA,decA);
    new_val = F(initA);
    if sfo_opt_get(opt,'verbosity_level',0)>0
        fprintf('Value after downward pass: \t%f\n',new_val);
    end
    if new_val<=cur_val+TOL
        break
    end
end
initAc = sfo_setdiff_fast(V,initA);
if F(initAc)>new_val
    if sfo_opt_get(opt,'verbosity_level',0)>0
        fprintf('Complement achieves higher value!\n');
    end    
    sset = initAc;
else
    sset = initA;
end

