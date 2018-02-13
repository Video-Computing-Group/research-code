% Andreas Krause (krausea@gmail.com)
% pSPIEL helper function: Compute cost and edges of a placement by (approximately) 
% solving steiner tree problem using the MST heuristic
%
% function  [cost,edges,steinernodes] = sfo_pspiel_get_cost(A,D,dists)
% A: set of nodes (indices in D)
% D: adjacency matrix
% dists: all pairs shortest path closure of D (optional, more efficient)
% cost: cost of the MST in dists connecting A. Is 2-approximation to
% minimum Steiner tree cost
% edges: edges in steiner tree
% steinernodes: additional nodes selected by expanding steiner tree
%
% Example: See tutorial script.

function  [cost,edges,steinernodes] = sfo_pspiel_get_cost(A,D,dists)
if ~exist('dists','var')
    dists = sfo_pspiel_sp(D);
end
% connect all selected nodes to a tree:
% This algorithm computes the MST on the selected nodes w.r.t. the SP
% metric; then expands edges by shortest paths. This is a (constant
% factor 2) approximation to the Steiner tree problem.
edges = [];
steinernodes = [];
nedges = 0;
% get MST on shortest path matrix
[cost,mstedges] = sfo_pspiel_mst(dists(A,A));
if size(mstedges,2)==2
    % non-trivial tree
    for j = 1:size(mstedges,1)
        % expand each MST edge in terms of the shortest paths
        [tmp, path] = sfo_pspiel_dijkstra(D, A(mstedges(j,1)), A(mstedges(j,2)));
        if length(path)>2 
            steinernodes = [steinernodes, path(2:(length(path)-1))];
        end
        for i=2:length(path)
            nedges = nedges+1;
            edges(nedges,:)=[path(i-1),path(i)];
        end
    end
    % make sure steinernodes are unique
    steinernodes = sfo_setdiff_fast(steinernodes,A);
else
    steinernodes = [];
end
