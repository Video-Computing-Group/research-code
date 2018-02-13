% Helper function for sfo_pspiel, by Andreas Krause (krausea@gmail.com)
% Example: See sfo_tutorial.m
% solve pSPIEL for fixed value of R
% last parameter pdallok controls whether it's ok to use the trivial PD or not
function [A, E, result] = sfo_pspiel_fixed_r(F,V,Q,D,R,dists,pdallok,Vroot)

result.failed = 0;

% selecting a padded decomposition
done = 0; trial = 0; maxTrials = 10;
while ~done && trial<maxTrials;
    trial = trial+1;
    [pd,done,a,successrate] = sfo_pspiel_pd(V,dists,R,1/2);
    if ~pdallok && length(pd)==1 && isempty(sfo_setdiff_fast(V,pd{1}))
        % we got a trivial partition. Let's avoid that by reducing R
        done = 0;
        R = R/2;
    elseif ~done
        % we didn't find a PD that achieves success probability >= 1/2.
        % Increase R.
        R = 2*R;
    end
end

% construct modular approximation graph (MAG)
[reward,adj,roots,gs,gimp,node2cluster] = sfo_pspiel_construct_mag(F,pd,dists,Vroot);

% for rooted k-MST, we build the MAG on the residual submodular function
% hence the reward of the root node is ignored --> need to reduce quota
if Vroot>0
    valRoot = F(Vroot);
else
    valRoot = 0;
end
% compute the Quota-MST solution on the MAG
allels = sfo_pspiel_kmst(reward,adj,Q-valRoot,roots);

if isempty(allels) 
    % failed to satisfy quota!
    A = []; E = [];
    result = prepareResult(1,[],[],[],0,1e99,0,Q,R,pd,a,successrate,gs,gimp);
    return    
else
    % transfering solution from MAG back to original graph
    selectednodes = sfo_pspiel_transfer_back(allels,gs,node2cluster);
    
    % now greedily augment to satisfy quota
    selectednodes = sfo_pspiel_augment(F,V,selectednodes,Q,dists);
    
    % solve Steiner tree problem to connect selectednodes
    [totalCost,selectededges,supportnodes] = sfo_pspiel_get_cost(selectednodes,D,dists);
    
    A = [selectednodes supportnodes];
    E = selectededges;
    
    selectedVal = F(selectednodes); totalVal = F(A);
        
    % generate result structure
    result = prepareResult(0,selectednodes,supportnodes,selectededges,totalVal,totalCost,selectedVal,Q,R,pd,a,successrate,gs,gimp);
end


%% Now we set up the modular approximation graph
function [reward,adj,roots,gs,gimp,node2cluster] = sfo_pspiel_construct_mag(F,pd,dists,Vroot)
nnodes = 0; nedges = 0; %these count the #nodes and #edges used

reward = []; % modular reward per node (takes indices in MAG)
node2cluster = []; %first column is cluster id (in pd), second column is node # in cluster acc. to greedy

edge = []; % indices from to in MAG
edgeweight = []; % edgeweights in MAG
nodeindex =  {}; % nodeindex{i} contains MAG indicies of nodes in cluster i
roots = []; % set of cluster centers, i.e., roots(i) = nodeindex{i}(1)
gs = {}; %greedy sequences (indices in V)
gimp = {}; %greedy improvements 

if Vroot>0
    % need to make sure root is part of MAG
    % first find cluster that contains root
    rootPDind = [];
    for i = 1:length(pd)
        if sum(pd{i}==Vroot)>0 && length(pd{i})>1
            rootPDind = i;
            break;
        end
    end
    % make root a separate cluster
    if ~isempty(rootPDind)
        % remove
        pd{rootPDind} = sfo_setdiff_fast(pd{rootPDind},Vroot);
    end
    pd{end+1} = Vroot;
    % now define a new residual submodular function F'(A) = F(A u root)-F(root)
    F = sfo_fn_residual(F,Vroot);
end

for i = 1:length(pd);
    % compute greedy subsets and modular approximation per cluster
    ni = length(pd{i}); %number of clusters
        
    if pd{i}(1) ~= Vroot || length(pd{i})>1
        [gs{i},scores] = sfo_greedy_lazy(F,pd{i},ni);
    else
        % in case we have a root, that will give value 0 and greedy will
        % not pick it.
        gs{i} = Vroot;
        scores = 1e-99;
    end
    if isempty(gs{i}) % require that the greedy algorithm picks at least 1 element per cluster. Maybe F == 0?
        error('greedy algorithm did not pick any elements!');
    end

    % compute greedy improvements F(Aj)-F(A_{j-1}) for i-th cluster
    gimp{i} = scores(1:end)-[0 scores(1:(end-1))];
    
    % compute increase in MST costs per cluster
    mstweight = zeros(1,ni);
    for j = 1:length(gs{i})
        mstweight(j) = sfo_pspiel_mst(dists(gs{i}(1:j),gs{i}(1:j)));

        % add a new node to MAG
        nnodes = nnodes+1;
        reward(nnodes) = gimp{i}(j); % assign greedy improvment score
        node2cluster(nnodes,:)=[i,j]; % i.e., MAG node nnodes is j-th node in i-th cluster
        
        if j>1 
            % not cluster center; need to attach node to center by an edge
            nedges = nedges+1;
            edge(nedges,:)=[nnodes-1,nnodes];
            
            % edge cost is difference in MST costs
            edgeweight(nedges) = max(0,mstweight(j)-mstweight(j-1)); 
            nodeindex{i} = [nodeindex{i},nnodes];
        else
            % new cluster center. Will be a possible root for the k-MST
            % algorithm
            nodeindex{i} = nnodes;
            roots = [roots,nnodes];
        end
        
    end
end
if Vroot>0
    roots = roots(end);
end

% now create the "core" of G', which is a complete subgraph (clique) connecting the
% cluster centers.
for i = 1:length(pd)
    for j = (i+1):length(pd)
        nedges = nedges+1;
        % connect centers of clusters i and j
        edge(nedges,:) = [nodeindex{i}(1),nodeindex{j}(1)];

        % edge weight is shortest path distance between cluster centers
        edgeweight(nedges) = dists(gs{i}(1),gs{j}(1)); 
    end
end

% generate adjacency matrix representation for k-MST algorithm
INFINITY = 1e99; %guardian value for non-existent edges
adj = INFINITY*ones(nnodes,nnodes);
for j=1:nedges
    adj(edge(j,1),edge(j,2))=edgeweight(j);
    adj(edge(j,2),edge(j,1))=edgeweight(j);
end
for i = 1:nnodes
    adj(i,i)=0;
end

%% Now transfer back from MAG to original graph
% A is the MAG nodes selected
function selectednodes = sfo_pspiel_transfer_back(A,gs,node2cluster)
% now we need to figure out how many nodes are selected in each cluster
numpercl = zeros(1,length(gs));
for y = A
    % for indexing
    clind = node2cluster(y,1); % the cluster index
    clpos = node2cluster(y,2); % position in the cluster
    numpercl(clind)=max(numpercl(clind),clpos);
end
selectednodes = []; % the real indices selected in V

% now select nodes in V: for each cluster i pick the first numpercl(i) greedy elements 
for i=1:length(gs)
    selectednodes = [selectednodes, gs{i}(1:numpercl(i))];
end


%% Greedily augment the network to achieve quota F(A)>=Q
% Due to (R,gamma)-locality, the selected nodes might achieve quota Q on
% the MAG, but less than Q on F. Here, we simply greedily add nodes until
% the quota is met.
function A = sfo_pspiel_augment(F,V,A,Q,dists)
while 1
    currentScore = F(A);
    if (currentScore>=Q) % satisfied quota
        break
    end
    Vc = sfo_setdiff_fast(V,A); % elements that are still available
    Fcur = init(F,A); % optimization for speed
    scores  = zeros(1,length(Vc)); % scores for adding each available element
    for i = 1:length(Vc)
        Fcur = init(Fcur,A);
        newCost = min(dists(A,Vc(i))); % cost of adding = cheapest shortest path connection cost
        newScore = max(Q,inc(Fcur,A,Vc(i))); % benefit of adding; don't care if > Q
        scores(i) = (newScore-currentScore)/max(newCost,1e-10); % score = benefit / cost ratio
    end
    [tmp best] = max(scores); % find best element
    A = [A Vc(best)];
end


%% prepare the result output structure
function result = prepareResult(failed,selectednodes,supportnodes,selectededges,totalVal,totalWeight,selectedVal,quota,R,pd,a,successrate,gs,gimp) %,reward,edge,edgeweight,adj)
result.failed = failed;
result.selectedNodes = selectednodes;
result.supportNodes = supportnodes;
result.selectedEdges = selectededges;
result.totalVal = totalVal;
result.totalWeight = totalWeight;
result.selectedVal = selectedVal;
result.parameters.quotaMI = quota;
result.parameters.localityR = R;
result.paddedDecomposition.clusters = pd;
result.paddedDecomposition.alpha = a;
result.paddedDecomposition.successRate = successrate;
result.greedySets = gs;
result.greedyImprovements = gimp;
