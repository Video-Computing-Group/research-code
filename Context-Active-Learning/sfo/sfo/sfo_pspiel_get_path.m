% Andreas Krause (krausea@gmail.com)
% pSPIEL Helper function: Convert tree into tour using the spanning
% tree heuristics (2 approx). Allows visiting the same node twice
% Based on implementation by Amarjeet Singh
%
% function [tour, E, cost] = sfo_pspiel_get_path(A,D) 
% A: set of nodes
% D: distance matrix (A indices subsets of rows and columns of D)
% Vroot: if Vroot is specified (which needs to be part of A), it
%   will be made first and last node in tour.
% tour: a sequence of nodes
% E: edge representation of tour
% cost: cost of tour
%
% returns a tour that's at most twice as long as the best tour
%
% Example: See tutorial script.

function [tour E cost]= sfo_pspiel_get_path(set,D,Vroot) 

if length(set)<2 %trivial solution
    tour = set;
    E = [];
    cost = 0;
    return
end

d1 = D(set,set);

% first get a minimum spanning tree
[tmp,T] = sfo_pspiel_mst(d1);
A = zeros(size(d1));
for i = 1:size(T,1)
    A(T(i,1),T(i,2))=1;
    A(T(i,2),T(i,1))=1;
end

order = sfo_pspiel_get_path_fast_dfs(A,1,[]);

% % identify nodes of odd degree
% oddNodes = find(mod(sum(A),2)>0);
% 
% % find a perfect matching of the odd degree nodes
% M1 = sfo_pspiel_get_path_perfect_matching(d1(oddNodes,oddNodes));
% 
% % combine MST and matching to get an Eulerian graph
% A(oddNodes,oddNodes) = A(oddNodes,oddNodes) + M1;
% 
% % get an eulerian tour 
% order = sfo_pspiel_get_path_euler_tour(A,1);

% get back original indices
tour = set(order);

% apply shortcutting
tour = sfo_pspiel_get_path_2opt(tour,D);

if exist('Vroot','var')
    if Vroot>0
        tour = tour(1:(end-1));
        rootpos = find(tour==Vroot);
        tour = tour([rootpos:end 1:rootpos]);
    end
end

% create edge representation
E = [tour(1:(end-1))', tour(2:end)'];
cost = 0;
for i = 1:size(E,1)
    cost = cost+D(E(i,1),E(i,2));
end

%% Do a DFS on a tree
function result = sfo_pspiel_get_path_fast_dfs(A,root,resid)
result = root;
resid = [resid root];
ch = sfo_setdiff_fast(find(A(root,:)),resid);
for i = 1:length(ch)
    subtour = sfo_pspiel_get_path_fast_dfs(A,ch(i),[resid i]);
    result = [result subtour root];
end


%% Find a minimum perfect matching on graph G with an even number of nodes
% Uses the matlab bintprog solver (often very fast). Could be replaced by the efficient
% deterministic algorithm
function [M,fval] = sfo_pspiel_get_path_perfect_matching(G)
%M is a binary matrix, where the elements with value 1 represent the matchings
%fval is the cost of the matching

n=max(size(G));

delta = zeros(n);
d = [];
c = 0;
for i=1:n-1
    for j=i+1:n
        c = c+1;
        d(c) = G(i,j);
        delta(i,j) = c;
        delta(j,i) = c;
    end
end

A = zeros(n,c);

for i=1:n
    for j=1:n
        if(i~=j)
            A(i,delta(i,j)) = 1;
        end
    end
end

[X,fval] = bintprog(d,zeros(1,c),0,A,ones(n,1),[],optimset('Display','off'));

M = zeros(n);

for i=1:n-1
    for j=i+1:n
        if(X(delta(i,j)) > 0)
            M(i,j) = 1;
            M(j,i) = 1;
        end
    end
end

%% Returns a vector representing an eulerian tour on graph D
function tour = sfo_pspiel_get_path_euler_tour(D, startNode)
%Input graph D is an integer valued matrix. Element D(i,j) indicates the
%number of distinct edges from node i to node j.
%D should be undirected and all nodes should be of even degree.
%startNode is the index of the node we want to start from the "to visit" set

tour=[];
currentNode = startNode;
flag = 1;
allDone = 0;

while(flag > 0)
    tour = [tour currentNode];
    row = D(currentNode,:);
    outgoing = find(row);
    
    if(min(size(outgoing)) < 1)
        error('Something went wrong: Probably some odd-degree vertex');
    else
        next = outgoing(1);
    
        D(currentNode,next) = D(currentNode,next) - 1;
        D(next,currentNode) = D(next,currentNode) - 1;
        
        currentNode = next;
        
        if(next == startNode)
            flag = 0;
            tour = [tour startNode];
        end
    end
end


while(allDone == 0)
    [I,J]=find(D);

    if(min(size(I)) < 1) 
        allDone = 1;
        %return; % this means we have used all edges
    else

    % if there are still unused edges...
    
    flag = 1;
    c = 1;
    
    while(flag > 0)
        currentNode = I(c);
        x = find(tour == currentNode);
        if(min(size(x)) < 1)
            c = c + 1;
            if(c > max(size(I)))
                error('Something went wrong: Probably disconnected graph');
            end
        else
            flag = 0;
        end
    end
    
    tourIndex = x(1);
    
    subtour = [];
    startNode = currentNode;
    flag = 1;
    
    while(flag > 0)
        subtour = [subtour currentNode];
        row = D(currentNode,:);
        outgoing = find(row);
        
        if(min(size(outgoing)) < 1)
            error('Something went wrong: Probably some odd-degree vertex');
        else
            next = outgoing(1);
            
            D(currentNode,next) = D(currentNode,next) - 1;
            D(next,currentNode) = D(next,currentNode) - 1;
            
            currentNode = next;
            
            if(next == startNode)
                flag = 0;
                subtour = [subtour startNode];
            end
        end
    end
    toursize = length(tour);
    subtoursize = length(subtour);
    tour = [tour(1:tourIndex) subtour(2:subtoursize) tour(tourIndex+1:toursize)];
    
    end
    
end

n = length(tour);
t=tour(1);


for i=2:n
    if(ismember(tour(i),tour(1:i-1)) < 1)
        t = [t tour(i)];
    end
end


tour = [t t(1)];

%% Shortcut the tour by applying 2-opt heuristic
function tour = sfo_pspiel_get_path_2opt(tour,D)
moveDone = 1;
while (moveDone == 1)
    moveDone = 0;
    for i = 1:(length(tour)-3)
        if(moveDone == 1)
            break;
        end
        for j = i+2:(length(tour)-1)
            origLength = D(tour(i),tour(i+1)) + D(tour(j),tour(j+1));
            newLength = D(tour(i),tour(j)) + D(tour(i+1),tour(j+1));
            if(newLength < origLength)
                moveDone = 1;
                tmpPath = tour;
                optIndex = [1:i,j:-1:(i+1),(j+1):length(tour)];
                tour = tmpPath(optIndex);
            end
        end
    end
end
i = 1;
while i<length(tour)
    if tour(i)==tour(i+1)
        tour = tour([1:i (i+2):end]);
    else
        i = i+1;
    end
end
