% Andreas Krause (krausea@gmail.com)
% Trade off utility  F and (tree) connection cost using 
% the pSPIEL algorithm of Krause et al. (IPSN '06)
%
% function [A, E, result] = sfo_pspiel(F,V,Q,D,opt)
% F: monotonic submodular function
% V: index set
% Q: quota
% D: Cost matrix. D(i,j) is edge cost from V(i) to V(j)
% opt (optional): parameter struct, depending on
%
% pspiel_max_iter: number of restarts (optional, default = 10)
% pspiel_root: if positive indicates root of tree. If <0, no root is chosen.
% pspiel_fixed_r: locality constant (F is assumed to be (R,gamma)-local (optional)
%
% returns set of nodes A and edges E as in IPSN '06, as well as additional
% information in result.
%
% Example: See tutorial script.


%%
function [A E result] = sfo_pspiel(F,V,Q,D,opt) 
if ~exist('opt','var')
    opt = sfo_opt;
end
maxIter = sfo_opt_get(opt,'pspiel_max_iter',10);
Vroot = sfo_opt_get(opt,'pspiel_root',-1);

% compute and cache shortest path metric
dists = sfo_pspiel_sp(D);

if isfield(opt,'pspiel_fixed_r')
    R = sfo_opt_get(opt,'pspiel_fixed_r');
    % run only a single trial with fixed R
    [A, E, result] = sfo_pspiel_fixed_r(F,V,Q,D,R,dists,0,Vroot);
else
    % we need to find a good value for the locality constant R
    numRs = 20;
    Rs = sfo_pspiel_get_r_range(dists,numRs);

    % best solution found so far
    bestCost = inf;
    
    % store solutions of random trials
    A = {}; E = {}; result = {}; val = zeros(1,maxIter); cost = zeros(1,maxIter); 
    
    % do multiple random trials
    for i = 1:maxIter
        if i == 1
            % in the first iteration, always pick maximum R (last in Rs)
            % This returns a padded decomposition with a single cluster
            % effectively doing Greedy Connect
            Rind = numRs;
        else
            % pick one of the values
            Rind = floor(rand*(numRs-1))+1;
        end
        R = Rs(Rind);
        
        % run pSPIEL for fixed value of R
        [A{i} E{i} result{i}] = sfo_pspiel_fixed_r(F,V,Q,D,R,dists, i==1,Vroot);
        
        % store result
        val(i) = result{i}.totalVal;
        cost(i) = result{i}.totalWeight;
        if val(i)>=Q && cost(i)<bestCost
            bestCost = cost(i);
        end
        if sfo_opt_get(opt,'verbosity_level',1)>0
           disp(sprintf('iteration %d: best cost = %f, R = %f, v = %f, c = %f',i,bestCost,R,val(i),cost(i)));
        end
    end
    if sum(val(:)>=Q)>0
        
        % found a feasible solution satisfying quota. Return result
        best = find(cost==min(cost(val(:)>=Q)),1);
        A = A{best}; E = E{best};  result = result{best};
    else
        % did not find feasible solution. Fail.
        if sfo_opt_get(opt,'verbosity_level',1)>0
            disp('Quota could not be satisfied!');
        end
        A = []; E = []; result = result{1,1};
    end
    
end
