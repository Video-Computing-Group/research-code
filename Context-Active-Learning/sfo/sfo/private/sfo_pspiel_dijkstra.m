% Andreas Krause (krausea@gmail.com)
% pSPIEL helper function: Compute shortest paths using dijkstra
% based on implementation by Xiaodong Wang
%
% function [distance, path, totalCost] = sfo_pspiel_dijkstra(adj, s, d)
% path: the list of nodes in the path from source to destination
% distances: all distances
% totalCost: cost of path
% adj: adjacency matrix
% s: source node index;
% d: destination node index;
%
% Example: See tutorial script.

function [distance, path, totalCost] = sfo_pspiel_dijkstra(adj, s, d)

n = size(adj,1);

farthestPreviousHop = zeros(1,n);
farthestNextHop = zeros(1,n);
for i = 1:n
    % initialize the farthest node to be itself;
    farthestPreviousHop(i) = i; 
    farthestNextHop(i) = i;
end;
    
% all the nodes are un-visited;
visited(1:n) = 0;

distance(1:n) = inf;    % it stores the shortest distance between each node and the source node;
parent(1:n) = 0;

distance(s) = 0;
for i = 1:(n-1),
    temp = zeros(1,n);
    temp(visited~=0)=inf;
    temp(visited==0)=distance(~visited);
     [t, u] = min(temp);    % it starts from node with the shortest distance to the source;
     visited(u) = 1;       % mark it as visited;
     toupdate = find(adj(u, :)+distance(u) < distance);
     for v = toupdate
         if ( ( adj(u, v) + distance(u)) < distance(v) )
             distance(v) = distance(u) + adj(u, v);   % update the shortest distance when a shorter path is found;
             parent(v) = u;                                     % update its parent;
         end;             
     end
end;

path = [];
if parent(d) ~= 0   % if there is a path!
    t = d;
    path = [d];
    while t ~= s
        p = parent(t);
        path = [p path];
        
        if adj(t, farthestPreviousHop(t)) < adj(t, p)
            farthestPreviousHop(t) = p;
        end;
        if adj(p, farthestNextHop(p)) < adj(p, t)
            farthestNextHop(p) = t;
        end;

        t = p;      
    end;
end;

totalCost = distance(d);

return;