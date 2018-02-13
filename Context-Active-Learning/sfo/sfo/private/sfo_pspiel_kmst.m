% Andreas Krause (krausea@gmail.com)
% pSPIEL helper function: Approximately solve Quota-MST problem on MAG
% This algorithm computes a log^3 n approximation to the quota-MST for a
% graph defined by the adjacency matrix adj, and the rewards per node
% defined by reward
% This is the Multiple-Kruskal like algorithm from Blum et al's STOC 95 paper.
%
% function [nodes,cost,rew] = sfo_pspiel_kmst(reward,adj,quota,roots)
% reward: vector of rewards for each node of adj
% adj: adjacency matrix
% quota: want sum(reward(nodes))>=quota
% roots: indices to consider as possible roots
% cost: cost of quota-MST
% rew: reward attained at solution
%
% Example: See tutorial script.

function [nodes,cost,rew] = sfo_pspiel_kmst(reward,adj,quota,roots)
% get and cache shortest path closure for adjacency matrix
sps = sfo_pspiel_sp(adj);

cost = []; nodes = {}; rew = []; % for storing solutions
for i = 1:length(roots)
    % try out every root
    [nodes{i}, cost(i), rew(i)] = sfo_pspiel_kmst_find_l(reward,adj,quota,roots(i),sps);
end
% pick best solution
best = find(cost == min(cost),1);
nodes = nodes{best}; cost = cost(best); rew = rew(best);


%% solve rooted k-MST problem for fixed root
function [nodes,cost, rew] = sfo_pspiel_kmst_find_l(reward,adj,quota,root,sps)

if quota<= reward(root)
    % trivial solution: only pick root
    nodes = root;
    cost = 0;
    rew = reward(root);
    return
end

% We need to find a bounding radius L that controls how far the k-MST can
% extend from the root. Need to search for this value. We do this by
% performing a halfing search, starting from an easy upper bound.

% get distances of all nodes from root
[dists, I] = sort(sps(root,:));
% cumulative reward of those nodes
cumrew = cumsum(reward(I));

if (cumrew(end)>=quota)
    % get easy upper bound Lmax on L
    Lmax = sum(dists(1:find(cumrew>=quota,1)));
    
    % store all solutions for different values of L
    cost = []; rew = []; nodes = {};
    for i = 1:10
        % try out new value of L
        L = Lmax*.5^(i-1);
        
        % solve k-MST for fixed value of L
        comps =  sfo_pspiel_kmst_fixed_l(reward,adj,quota,root,L,sps);
        
        nodes{i} = sfo_unique_fast([comps{:}]); % nodes selected
        rew(i) = sum(reward(nodes{i}));
        if rew(i)<quota
            % radius L is too small to satisfy quota
            break
        end
        % store cost of solution
        cost(i) = sfo_pspiel_mst(sps(nodes{i},nodes{i}));
    end
    if length(cost)>0
        % found solution: return
        best = find(cost==min(cost),1);
        nodes = nodes{best};
        cost = cost(best);
        rew = rew(best);
        return
    end
end

% failed to satisfy quota
nodes = [];
cost = inf;
rew = 0;


%% solve k-MST problem using the Blum et al. algorithm for fixed root and fixed radius L
function comps = sfo_pspiel_kmst_fixed_l(reward,adj,quota,root,L,sps)

dbg = false; % debug output
k = quota;
reward = max(reward,1e-10); % guarantee that the reward is positive
n = length(reward); 

% will only be allowed to pick non-marked vertices
marked = zeros(1,n);
marked(sps(root,:)>L)=1; % mark all nodes further than L from root

allcomp=[];
iter = 0;
comps = {};
% This loop is called at most log k times
while k>0 && sum(marked==1)<n
    iter = iter + 1;
    % call the bicriterion optimization algorithm to generate a cluster
    comp = sfo_pspiel_kmst_merge_cluster(reward,marked,adj,k/4,sps);
    
    % mark nodes in the returned cluster
    marked(comp)=1;
    
    % combine all found components
    allcomp = sfo_unique_fast([allcomp,comp]);
    
    % figure out how much quota we still need to satisfy
%     if sum(allcomp==root)==0
%         quotatogo = quota;
%     else
%         quotatogo = quota-reward(root);
%     end
    k = quota - sum(reward(sfo_unique_fast([root allcomp])));
    
%    k = quotatogo - sum(reward(allcomp));
    if dbg
        disp(sprintf('Found new cluster, quota to go: %f',k));
    end
    comps{iter}=comp;
end
if sum(allcomp==root)==0 % if we didn't pick the root, add it
    comps{end+1} = root;
end

%%
% This procedure, which is part of Blum et al's algorithm, is Kruskal-like. 
% It grows clusters until one is large enough to reach the quota. 
% It will only consider unmarked nodes. The distance between clusters 
% is the cost/benefit ratio: edgecost/(smaller weight of both clusters)
function comp=sfo_pspiel_kmst_merge_cluster(reward,marked,adj,quota,sps)
n = length(reward);
infty = 1e99;

comps = {}; % contains all components
for j = 1:n
    % initialize components with a single element
    comps{j}.els=[j];
    comps{j}.active = 1;
end

d=infty*ones(n,n); % distances between clusters
maxVal=-infty;
argmaxVal = -1;
markedSet = find(marked ==1); % indices of marked nodes

for c1 = 1:n % iterate through the components
    if marked(c1) %ignore marked vertices
        comps{c1}.active = 0;
        continue
    end
    % get reward
    mi = reward(c1);
    comps{c1}.reward = mi;
    if mi>maxVal % best reward so far
        maxVal = mi;
        argmaxVal = c1;
    end

    dists = sps(c1,:); %shortest path distances to everyone else
    for c2 = (c1+1):n
        if marked(c2)
            continue
        end
        dist = dists(c2);
        marg = min(reward(c1),reward(c2));
        d(c1,c2)=dist/marg; % initialize greedy rule matrix
    end
end  

%d is the greedy rule matrix
activeComps = n-sum(marked);
if maxVal>quota || activeComps==1
    % found a singleton solution, or have only one active component
    comp = argmaxVal;
    return;
end
while activeComps>1
   
    % find pair of clusters minimizing cost/benefit ratio
    [c1,c2] = ind2sub(size(d),find(d == min(d(:)),1));
    
    % compute the path connecting the components
    path = sfo_pspiel_kmst_find_connection(c1,c2,comps,adj,sps);
    
    % put all elements from c2 into c1. Also add path connection
    comps{c1}.els = sfo_unique_fast([comps{c1}.els,comps{c2}.els,path]);
    comps{c1}.reward = sum(reward(sfo_setdiff_fast( comps{c1}.els, markedSet) ));
    % clear and disable cluster c2
    comps{c2}.els = [];
    comps{c2}.active = 0;
    comps{c2}.reward = 0;
    activeComps = activeComps-1;
    d(:,c2)=infty;
    d(c2,:)=infty;
    
    % check if we're done
    if comps{c1}.reward>=quota
        break
    end
    
    % now update distances for greedy rule
    for i=[1:(c1-1) (c1+1):n]
        if comps{i}.active == 0
            continue
        end
        
        minGain = min(comps{c1}.reward,comps{i}.reward); % compare rewards
        minDist = min(min(sps(comps{c1}.els,comps{i}.els))); % get connection cost
        d(i,c1) = minDist/minGain;
    end
end
comp = comps{c1}.els;

%% helper function that returns the shortest path and distance to connect two clusters.
function path = sfo_pspiel_kmst_find_connection(c1,c2,comps,adj,sps)
v1 = comps{c1}.els; 
v2 = comps{c2}.els;

% Get the two vertices s1 and s2 that minimize the cluster distance
spsSel = sps(v1,v2);
[s1 s2] = ind2sub(size(spsSel),find(spsSel==min(spsSel(:)),1));

% use the shortest path matrix to recover the shortest path
path = sfo_pspiel_kmst_recoverSP(s1,s2,adj,sps);

% only return middle nodes of the path
path = path(2:(length(path)-1));

%% find shortest path betweet s1 and s2, using shortest path matrix
function path = sfo_pspiel_kmst_recoverSP(s1,s2,adj,sps)
path = s1;
s = s1;
while s~=s2 % not arrived yet
    % dynamic programming trick to find next node on path
    ndist = adj(s,:)+sps(:,s2)'; ndist(path)=inf;
    % get best next node and add to path
    s = find(ndist == min(ndist),1);
    path = [path s];
end
