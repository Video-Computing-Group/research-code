% Implements Greedy splitting by Zhao et al
% Repeatedly uses Queyranne's algorithm to find the optimal two-partition
% minimizing the energy E(A)+E(V\A)
% author: Andreas Krause (krausea@gmail.com)
%
% function P = sfo_greedy_splitting(E,V,k)
% E is the energy function per cluster
% V is the index set
% k is the number of clusters
% 
% Returns partition P; factor 2 approximation to minimum energy k-partition
%
% Example: see tutorial script sfo_tutorial.m

function P = sfo_greedy_splitting(E,V,k)

P = {V};
for i = 1:k-1
    cand = {}; scores = inf*ones(1,length(P));
    % now iterate through all clusters that contain at least 2 elements
    for j = 1:length(P)
        disp(sprintf('iteration %d, cluster %d',i,j));
        V1 = P{j}; 
        if length(V1)==1
            continue;
        end
        Etotal = E(V1);
        % now symmetrize the energy function
        Fsym = sfo_fn_wrapper(@(A) E(A)+E(sfo_setdiff_fast(V1,A))-Etotal);
        % Find best split using Queyranne's algorithm
        A1 = sfo_queyranne(Fsym,V1);
        cand{j} = {A1, sfo_setdiff_fast(V1,A1)};
        scores(j) = Fsym(A1);
    end
    % now greedily pick the best candidate partition
    argmin = find(scores==min(scores),1);
    P = [{P{1:(argmin-1)}} cand{argmin} {P{(argmin+1):length(P)}}];
end
