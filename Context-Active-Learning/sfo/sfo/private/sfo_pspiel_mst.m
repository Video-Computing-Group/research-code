% Andreas Krause (krausea@gmail.com)
% pSPIEL helper function: Compute a minimum spanning tree (MST)
% based on the implementation of F. van den Berg
%
% function [weight,Xmst] = sfo_pspiel_mst(D)
% D: adjacency matrix
% weight: cost of MST connecting all nodes in D
% Xmst: edges of the MST
%
% Example: See tutorial script.

function [weight,Xmst] = sfo_pspiel_mst(D)

nX = size(D,1);
Dmax = max(max(triu(D,1)))*10;

[Dmin,Dwin] = min(D(1,2:nX));

Xmst(1,:) = [1 Dwin+1];
if isempty(Dmin)
    weight = 0;
else
    weight = Dmin;
end
for a=2:nX-1
    mindist = Dmax;
    Xmstlist = sfo_unique_fast(Xmst(:));
    Xmstnotlist = sfo_setdiff_fast(1:nX,Xmstlist);
    for aa=Xmstlist
        [Dmin,Dwin] = min(D(aa,Xmstnotlist));
        if (Dmin < mindist)
            minindex = [aa Xmstnotlist(Dwin)];
            mindist = Dmin;
        end
    end
    Xmst(a,:) = minindex;
    weight = weight + mindist;
end
