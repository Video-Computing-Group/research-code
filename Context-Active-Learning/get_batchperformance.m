function [labeleddata, acc, n_label, newedgepot, newmodel] = get_batchperformance(labeleddata, batchdata, testset, link, edgepot, model)

[~,~, pmf] = svmpredict(ones(length(batchdata(:,end)),1), batchdata(:,2:end-1), model, '-b 1 -q');
prob = [];
for i = 1:length(model.Label)
    prob(:,i) = pmf(:,model.Label==i);
end

% Create graph
nodepot = prob;
n_nodes = size(batchdata,1);
adj = zeros(n_nodes);
W = [];
for i = 1:n_nodes
    dist = repmat(batchdata(i,2:end-1),[size(batchdata,1) 1])-batchdata(:,2:end-1);
    dist = sum(dist.^2, 2);
    dist(i) = Inf;
    W = [W; dist'];
    [~,pos] = sort(dist);
    pos = pos(1:5);
    adj(i,pos) = 1;
    adj(pos,i) = 1;
end
[maxP,pos] = max(prob);
indWeak = [find(maxP>0.9) pos(maxP>0.9)];

% Infer on the graph
edgestruct = UGM_makeEdgeStruct(adj, model.nr_class);
edgepot_stack = repmat(edgepot + 100, [1,1,edgestruct.nEdges]);
[nodebel, edgebel] = UGM_Infer_LBP(nodepot, edgepot_stack, edgestruct);

% Compute graph stats
H = get_entropy(nodebel');
M = get_mutualinformation(nodebel',edgebel, n_nodes, edgestruct.edgeEnds);

% Optimize to obtain the informative nodes for manual labeling
M = M + M';
SFMinData.Q = -M;
SFMinData.f = sum(M,2) - H';
SFMinData.gamma = 0;
save SFMinData SFMinData  % This is loaded in the SFO toolbox code

A = sfo_min_norm_point(sfo_fn_subsetSel,1:size(M,1));
F = sfo_fn_subsetSel;
SFMinData.gamma = 1.1*abs(F(A)/size(M,1));
save SFMinData SFMinData
ind = sfo_min_norm_point(sfo_fn_subsetSel,1:size(M,1));
ind = setdiff(ind,indWeak(:,1));

% Update models
labeleddata = [labeleddata; batchdata(ind,:);[batchdata(indWeak(:,1),1:end-1) indWeak(:,2)]];
newmodel = svmtrain(labeleddata(:,end), labeleddata(:,2:end-1), '-t 0 -c 1 -b 1 -q');
newedgepot = get_edgepot(edgepot, link, batchdata(ind,1), batchdata(ind,end));
newedgepot = get_edgepot(newedgepot, link, batchdata(indWeak(:,1),1), indWeak(:,2));

% Obtain batch performance
[~, acc, ~] = svmpredict(testset(:,end), testset(:,2:end-1), newmodel, '-b 1 -q');
n_label = length(ind) + size(indWeak,1);
