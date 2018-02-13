% Andreas Krause (krausea@gmail.com)
% Algorithm for informative path planning, based on pSPIEL
% algorithm by Krause et al. (IPSN '06)
%
% function [A, E, result] = sfo_pspiel_orienteering(F,V,B,D,maxIter,R)
% F: monotonic submodular function
% V: index set
% B: budget
% D: Cost matrix. D(i,j) is edge cost from V(i) to V(j)
% maxIter: number of restarts (optional, default = 10)
% Vroot: if positive indicates starting point of tour. If <0, no root is chosen.
%
% returns a sequence of nodes A, making a tour in the graph.
% pSPIEL will attempt to find a path maximizing F(A), subject
% to the constraint that the path cost is bounded by the budget B.
%
% Example: See tutorial script.


%%
function [A E result] = sfo_pspiel_orienteering(F,V,B,D,maxIter,Vroot) 
if ~exist('maxIter','var')
    maxIter = 20;
end
if ~exist('Vroot','var')
    Vroot = -1;
end

% compute and cache shortest path metric
dists = sfo_pspiel_sp(D);


% we need to find a good value for the locality constant R
numRs = 20;
Rs = sfo_pspiel_get_r_range(dists,numRs);

maxVal = F(V); % upper bound on best value;

% best solution found so far
bestVal = 0;

% store solutions of random trials
A = {}; E = {}; result = {}; val = zeros(1,maxIter); cost = zeros(1,maxIter); 

% do multiple random trials
testQ = maxVal/2; %bestVal+rand(1)*(maxVal-bestVal); %guess quota
fact = 2;
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
    [tree edges result{i}] = sfo_pspiel_fixed_r(F,V,testQ,D,R,dists, i==1,Vroot);
    [A{i} E{i} cost(i)]= sfo_pspiel_get_path(tree,D,Vroot);

    % store result
    val(i) = F(sfo_unique_fast( A{i})); 
    if val(i)>bestVal && cost(i)<=B
        bestVal = val(i);
        testQ = bestVal*1.05;
        fact = 1;
    else
        fact = fact*2;
        testQ = bestVal+(maxVal-bestVal)/fact; %guess quota
    end
    
    disp(sprintf('iteration %d: testq = %f, best value = %f, R = %f, v = %f, c = %f',i,testQ,bestVal,R,val(i),cost(i)));
    
end
if sum(cost(:)<=B)>0

    % found a feasible solution satisfying quota. Return result
    best = find(val==max(val(cost(:)<=B)),1);
    A = A{best}; E = E{best};  result = result{best};
else
    % did not find feasible solution. Fail.
    disp('Budget could not be satisfied!');
    A = []; E = []; result = result{1,1};
end
    
