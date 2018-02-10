function [consolidated_el] = remove_edges_from_consolidated_list(consolidated_el, affected_edges)

t1 = consolidated_el(1,1);
t2 = consolidated_el(1,2);
z1 = consolidated_el(1,3);
z2 = consolidated_el(1,4);
n1 = consolidated_el(1,5);
n2 = consolidated_el(1,6);
consolidated_el(1, :) = [];
rownum = find(consolidated_el(:,1)==t1&consolidated_el(:,2)==t2&consolidated_el(:,3)==z1&consolidated_el(:,4)==z2 & (consolidated_el(:,5) == n1 | consolidated_el(:,6) == n2));
if ~isempty(rownum)
    consolidated_el(rownum, :) = [];
end;
rownums = [];
if ~isempty(affected_edges)
    [~, rownums, ~] = intersect(consolidated_el(:,1:6), affected_edges, 'rows');
end;
if ~isempty(rownums)
    for r = rownums'
        t1 = consolidated_el(r,1);
        t2 = consolidated_el(r,2);
        z1 = consolidated_el(r,3);
        z2 = consolidated_el(r,4);
        n1 = consolidated_el(r,5);
        n2 = consolidated_el(r,6);
        consolidated_el(r, :) = [];
        rr = find(consolidated_el(:,1)==t1&consolidated_el(:,2)==t2&consolidated_el(:,3)==z1&consolidated_el(:,4)==z2 & (consolidated_el(:,5) == n1 | consolidated_el(:,6) == n2));
        if ~isempty(rr)
            consolidated_el(rr, :) = [];
        end;
    end;
end;

return;