% Andreas Krause (krausea@gmail.com)
% solve the submodular coverage problem using greedy algorithm,
% i.e., for additive cost function, finds (approximately) cheapest set that 
% achieves F(A)>=Q for some quota Q.
%
% function [A, stat] = sfo_cover(F,V,Q,opt)
% F: submodular function
% V: index set
% Q: quota minimum value to achieve
% opt (optional): option structure, depending on
%
% cost: Cost vector. i-th element is cost of element V(i)
% cover_max_cost: maximum cost allowed
% cover_tolerance: tolerance for stopping
%
% returns solution sset with F(sset)>=Q, or sset = V
% stat == 1 iff F(sset)>=Q, 0 otherwise
%
% Example: A = sfo_cover(sfo_fn_example,1:2,2);

function [A, stat] = sfo_cover(F,V,Q,opt)
if ~exist('opt','var')
    opt = sfo_opt;
end
% if no argument specified, assume unit cost
C = sfo_opt_get(opt,'cost',ones(1,length(V)));
maxCost = sfo_opt_get(opt,'cover_max_cost',inf);
TOL = sfo_opt_get(opt,'cover_tolerance',1e-6);
opt.greedy_cost_benefit = 1.0;

Fc = sfo_fn_trunc(F,Q);

% use the greedy algorithm on the truncated function
A = sfo_greedy_lazy(Fc,V,maxCost,opt);
if F(A)<Q-TOL
    stat = 0;
else
    stat = 1;
end
