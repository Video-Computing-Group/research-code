function [adj] = compute_adjacency_matrix(cord)

%#codegen
coder.inline('None');

adj = zeros(length(cord), length(cord));

for c = 1:length(adj)
    nbrs = find_neighboring_cells(cord, c, 'circle');
    for nbr = nbrs
        adj(c, nbr) = 1;
        adj(nbr, c) = 1;
    end;
end;

return;