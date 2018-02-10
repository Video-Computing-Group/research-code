function [c] = update_corr_matrix(c, affected_edges)

for i = 1 : size(affected_edges, 1)
    t1 = affected_edges(i,1);
    t2 = affected_edges(i,2);
    z1 = affected_edges(i,3);
    z2 = affected_edges(i,4);
    n1 = affected_edges(i,5);
    n2 = affected_edges(i,6);
    c{t1,t2}{z1,z2}(n1,n2) = 1;
end;

return;