% Andreas Krause (krausea@gmail.com)
% pSPIEL helper function: Computes the all-pairs shortest path solution 
% from distance matrix D, by repeatedly calling Dijkstra's algorithm n times
%
% function result = sfo_pspiel_sp(D)
% D: adjacency matrix
% result: shortest path closure matrix
%
% Example: See tutorial script.

function result = sfo_pspiel_sp(D)
n = size(D,1);
result = zeros(n,n);
for i = 1:n
    result(i,:) = sfo_pspiel_dijkstra(D,i,i);
end

