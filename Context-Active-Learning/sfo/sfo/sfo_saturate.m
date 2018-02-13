% Andreas Krause (krausea@cs.cmu.edu)
% outer binary search loop for saturate algorithm
%
% function [ssetMax,ssetMin,maxc,minc] = sfo_saturate(F,V,k,task,opt)
% F: thresholdable submodular function (takes set and threshold)
% V: index set
% k: #elements to pick
% task: 'maxthresh' or 'minthresh' for maximizing threshold
% opt: optional parameters; depends on:
%
% saturate_maxbound: bound on maximum value
% saturate_coverfn: algorithm used for coverage
% saturate_maxiter: max # iterations
% binsearch_tolerance: termination tolerance
%
% SATURATE builds on thresholdable functions;
% Implemented examples include sfo_fn_lincomb (useful for maximizing
% the minimum over a collection of submodular functions) and 
% sfo_fn_varred (useful for minimizing the worst-case remaining variance)
%
% Example: (also see sfo_tutorial.m)
%
% Suppose we want to maximize the minimum of functions F1 and F2:
% F1 = sfo_fn_infogain(merced_data.sigma,1:86,.1)
% F2 = sfo_fn_mi(merced_data.sigma,1:86)
% F = sfo_fn_lincomb({F1,F2},[1 1])
% sfo_saturate(F,1:86,10,'maxthresh')
%
% Suppose we want to minimize the worst-case remaining variance:
% Here, we use 'minthresh' instead
% F = sfo_fn_varred(merced_data.sigma,1:86)
% sfo_saturate(F,1:86,10,'minthresh')

function [ssetMax,ssetMin,maxc,minc] = sfo_saturate(F,V,k,task,opt)
if ~exist('opt','var')
    opt = sfo_opt;
end
maxc = sfo_opt_get(opt,'saturate_maxbound',sfo_maxbound(F,V)); % default: compute unconstrained bound
if strcmp(task,'maxthresh')
    direction = -1;
elseif strcmp(task,'minthresh')
    direction = 1;
else
    error('SATURATE: task must be maxthresh or minthresh');
end
coverFn = sfo_opt_get(opt,'saturate_coverfn', @(F, V, B) sfo_greedy_lazy(F, V, B));
TOL = sfo_opt_get(opt,'binsearch_tolerance',1e-6)*maxc;
maxiter = sfo_opt_get(opt,'saturate_maxiter',50);

minc = 0;
ssetMin = [];
ssetMax = [];

iter = 0;
while ~((length(ssetMax)<=k && iter>=maxiter) || abs(maxc-minc)<TOL)
    iter = iter+1;
    testc = mean([minc,maxc]);
    if sfo_opt_get(opt,'verbosity_level')>0
        disp(sprintf('interval [%f,%f], size(ssetMax,Min)=%d,%d',minc,maxc,length(ssetMax),length(ssetMin)));
    end

    satF = trunc(F,testc); 
    sset = coverFn(satF,V,k+1);
    
    if direction*length(sset)>direction*k
        minc = testc;
        ssetMin = sset;
    else
        maxc = testc;
        ssetMax = sset;
    end
    if (iter>maxiter && direction==-1 && length(ssetMin)<=k)
        ssetMax = ssetMin;
        break;
    end
end
