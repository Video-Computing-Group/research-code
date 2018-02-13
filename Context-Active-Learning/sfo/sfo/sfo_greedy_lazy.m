% Andreas Krause (krausea@gmail.com)
% compute the greedy subset selection using lazy evaluations
%
% function [sset,scores,evalNum] = sfo_greedy_lazy(F,V,B,opt) 
% F: submodular function
% V: index set
% B: budget 
% opt: optional parameter struct, with following entries:
%
% cost: Cost vector (optional). i-th element is cost of element V(i)
% greedy_use_cost_benefit: use cost benefit rule to select element
% greedy_check_indep: function to check if a set is independent
%   for optimizing F subject to a matroid constraint
% greedy_initial_sset: start subset
%
% returns solution sset with C(sset)<=B
% scores(i) is the greedy score after element i was added
% evalNum(i) is the # of function evaluations in iteration i
%
% Example: A = sfo_greedy_lazy(@sfo_fn_example,1:2,1)

function [sset,scores,evalNum] = sfo_greedy_lazy(F,V,B,opt) 
if ~exist('opt','var')
    opt = sfo_opt;
end

n=length(V);
C = sfo_opt_get(opt,'cost',ones(1,n));
useCB = sfo_opt_get(opt,'greedy_use_cost_benefit',0);
checkIndep  = sfo_opt_get(opt,'greedy_check_indep',@(A) 1);

deltas = inf*ones(1,length(V)); %initialize optimistically

TOL = sfo_opt_get(opt,'greedy_tolerance',1e-6);

sset = sfo_opt_get(opt,'greedy_initial_sset',[]); %start with empty set or specified

curVal = 0;
curCost = 0; %no cost

evalNum = []; scores = []; %keep track of statistics
i = 0;
while 1
    i = i+1;
    if sfo_opt_get(opt,'verbosity_level',0)>1
        fprintf('Iteration: %d\n',i);
    end
    bestimprov = 0;
    evalNum(i) = 0;
    
    F = init(F,sset);
    
    deltas(curCost+C>B) = -inf; %cannot afford
    [tmp,order] = sort(deltas,'descend');
    
    % Now let's lazily update the improvements
    for test = order
        if ~checkIndep([sset V(test)])
            deltas(test) = -inf;
        end
        if deltas(test)>=bestimprov % test could be a potential best choice
            evalNum(i) = evalNum(i) + 1;
            improv = inc(F,sset,V(test))-curVal;
            if useCB
                improv = improv/C(test);
            end
            deltas(test)=improv;
            bestimprov = max(bestimprov,improv);
        elseif deltas(test)>-inf
            break;
        end
    end
    argmax = find(deltas==max(deltas),1); % find best delta
    if deltas(argmax)>TOL % nontrivial improvement by adding argmax
        F = init(F,[sset V(argmax)]);
        sset = [sset,V(argmax)];
        if (useCB) %need to account for cost-benefit ratio
            curVal = curVal + deltas(argmax)*C(argmax);
        else
            curVal = curVal + deltas(argmax);
        end
        curCost = curCost + C(argmax);
        scores(i) = curVal;
    else 
        break
    end
end
