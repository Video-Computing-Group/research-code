% Andreas Krause (krausea@gmail.com)
% pSPIEL helper function: Compute a padded decomposition.
% This implementation of an algorithm by A. Gupta  et al. (STOC '03) will return a 
% padded decomposition such that all clusters C in cl guarantee:
% 1) diam(cl)<a*R
% 2) every node in C is R padded with prob. at least successrate > padrate
% It does the necessary resampling etc.
%
% function [cl,done,a,successrate] = sfo_pspiel_pd(V,pdists,R,padrate)
% V: indices of nodes
% pdists: matrix encoding a metric satisfying triangle inequality
% R: padding radius
% padrate: minimum required successrate (default = 0.5)
%
% Example: See tutorial script.

function [cl,done,a,successrate] = sfo_pspiel_pd(V,pdists,R,padrate)
% select relevant subset V
dists = pdists(V,V);
if ~exist('padrate','var')
    padrate = 0.5;
end

n = size(dists,1);
maxiter = 4*log(n); % from STOC paper
done = false;

cl ={}; % initially, clusters are empty

for i = 1:maxiter
    a = 10^(rand(1)); %randomly pick a guess a in [1,10]
    
    % sample a padded decomposition with bound a
    cl = sfo_pspiel_pd_sample(dists,R,a);
    
    % filter the non padded nodes
    cl = sfo_pspiel_pd_filter(dists,R,cl);
    totlen = length([cl{:}]);
    successrate = totlen/n; % compute success rate
    if successrate>padrate
        done = true;
        break;
    end
end

for i = 1:length(cl)
    % map back to original indices
    cl{i} = V(cl{i});
end

%% sample a padded decomposition
% implements the algorithm of A. Gupta et al's STOC 2003 paper
% the ratio between cluster radius and likely paddedness is 64 dim(X)
% parameter a will override this scaling parameter
function clusters = sfo_pspiel_pd_sample(dists,R,a)
n = size(dists,1);
R = R*a; %required to make 4 R paddedness likely
R = R/4; %this is required to scale constants appropriately (see STOC paper)
% create an epsilon-Net with epsilon = R
rnet = sfo_pspiel_pd_create_rnet(dists,R);
rnet = rnet(randperm(length(rnet)));
rad = 2*R - rand(1)*R;
candidates = 1:n;
iter = 0;
clusters = {};
for y = rnet
    newcluster = candidates(dists(y,candidates)<=rad);
    if ~isempty(newcluster)
        iter = iter + 1;
        clusters{iter} = newcluster;
        candidates = setdiff(candidates,clusters{iter});
        if isempty(candidates)
            break
        end
    end
end

%% returns a new set of clusters such that all elements are R-padded
function result = sfo_pspiel_pd_filter(dists,R,pd)
n = size(dists,1);
allInds = 1:n;
iter = 0;
result = {};
for i = 1:length(pd)
    cl = [];
    for y = pd{i};
        ball = allInds(dists(y,:)<R);
        if isempty(sfo_setdiff_fast(ball,pd{i}))
            cl = [cl,y];
        end
    end
    if ~isempty(cl)
        iter = iter+1;
        result{iter}=cl;
    end
end

%% greedily creates R-Net (i.e. subset of nodes which are at least R apart but R-cover all other nodes);
function rnet = sfo_pspiel_pd_create_rnet(dists,R)
n = size(dists,1);
rnet = 1;
allInds = 1:n;
while true
    d2net = min(dists([1 rnet],:));
    notcovered = allInds(d2net>R);
    if isempty(notcovered)
        break;
    end
    rnet = [rnet, notcovered(1)];
end
