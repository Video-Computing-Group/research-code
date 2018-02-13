% The submodular-supermodular procedure of Narasimhan & Bilmes 
% Implemented by Andreas Krause (krausea@gmail.com)
% This algorithm is guaranteed to converge to a local optimum
%
% function A = sfo_sssp(F,G,V,opt)
% F: submodular function
% G: submodular function
% V: index set
% Returns a locally optimal solution to the problem A = argmin_A F(A)-G(A)
%
% Example: See sfo_tutorial.m

function A = sfo_ssp(F,G,V,opt)
if ~exist('opt','var')
    opt = sfo_opt;
end
TOL = sfo_opt_get(opt,'ssp_tolerance',1e-6);

N = length(V);
pi = V(randperm(N));
bestVal = inf;
A = [];
while 1
    Hw = sfo_ssp_modular_approx(G,pi);
    H = sfo_fn_wrapper(@(A) sum(Hw(sfo_unique_fast(A))));
    
    FM = sfo_fn_lincomb({F,H},[1,-1]);
    A = sfo_min_norm_point(FM,V,sfo_opt({'minnorm_init',A}));
    curVal = FM(A);
    D = sfo_setdiff_fast(V,A);
    D = D(randperm(length(D)));
    pi = [A D];
    if curVal<bestVal-TOL
        bestVal = curVal;
    else
        break;
    end    
end

function H = sfo_ssp_modular_approx(G,pi)
H = zeros(1,max(pi(:)));
W = [];
oldVal = G(W);
for i = 1:length(pi)
    newVal = G([W pi(i)]);
    H(pi(i)) = newVal-oldVal;
    oldVal = newVal;
    W = [W pi(i)];
end

