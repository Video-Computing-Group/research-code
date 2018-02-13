% The CELF algorithm from Leskovec et al, KDD '07
% Author: Andreas Krause (krausea@gmail.com)
% It performs cost benefit optimization, by simultaneously running the
% cost-benefit greedy and the cost-agnostic greedy algorithm;
% the better of the two solution is guaranteed to obtain a constant
% fraction .5*(1-1/e) approximation.
%
% function A = sfo_celf(F,V,B,opt)
% F: submodular function
% V: index set
% B: budget 
% opt (optional): parameter struct, depending on 
%
% cost: Cost vector. C(i) is cost of element V(i)
%
% Returns:
% A: Near optimal cost-benefit tradeoff solution
%
% Example: A = sfo_celf(@sfo_fn_example,1:2,1)

function A = sfo_celf(F,V,B,opt)
if ~exist('opt','var')
    opt = sfo_opt;
end
% run greedy algorithm ignoring costs
opt.greedy_use_cost_benefit = 0;
[A1,S1] = sfo_greedy_lazy(F,V,B,opt);
% run greedy cost-benefit greedy algorithm 
opt.greedy_use_cost_benefit = 1;
[A2,S2] = sfo_greedy_lazy(F,V,B,opt);
if sfo_opt_get(opt,'verbosity_level',0)>0
    disp(sprintf('Unit cost solution = %f, Cost-benefit solution = %f',S1(end),S2(end)));
end
% pick better of the two solutions
if S1(end)>S2(end)
    A = A1;
else
    A = A2;
end
