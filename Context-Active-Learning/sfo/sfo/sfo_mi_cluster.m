% Example program for clustering using mutual information
% Author: Andreas Krause (krausea@gmail.com)
%
% Greedily picks k "most informative" representatives
% Then finds the mutual information mincuts separating the cluster centers
% returns the partition P into clusters, and centers Agreedy
%
% function [P,Agreedy]=sfo_mi_cluster(sigma,k)
% sigma: covariance matrix
% x: numbder of clusters

function [P,Agreedy]=sfo_mi_cluster(sigma,k)
V = 1:size(sigma,1);
F = sfo_fn_mi(sigma,1:size(sigma,1));

% Agreedy contains representatives of each cluster
Agreedy = sfo_greedy_lazy(F,V,k);

P = {V};

for i = 1:(k-1)
    for j = 1:length(P)
        V1 = P{j};
        F = sfo_fn_mi(sigma,P{j});
        reps = intersect(V1,Agreedy);
        if length(reps)==1 % contains only a single cluster center
            continue
        end
        s = reps(1);
        t = reps(2);
        
        A = sfo_s_t_mincut(F,P{j},s,t);
        P = [P{1:(j-1)} {A, sfo_setdiff_fast(P{j},A)} P{(j+1):length(P)}];
        break;
    end
end
