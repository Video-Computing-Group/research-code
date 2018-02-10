function [edgelist] = compute_edge_list_from_adj(adj)

edgelist = [];

for i = 1:size(adj,1)
    J = find(adj(i,:)==1);
    edgelist = [edgelist; [i.*ones(length(J),1), J']];
end;

edgelist = sort(edgelist,2);
edgelist = unique(edgelist, 'rows');

return;
    