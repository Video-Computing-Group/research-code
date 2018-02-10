function [c] = data_association_atomic_quartet(el, c)

consolidated_el = [];
consolidated_el = [consolidated_el; [ones(size(el{1,1}{1,2},1),1)*[1 1 1 2], el{1,1}{1,2}]];
consolidated_el = [consolidated_el; [ones(size(el{2,2}{1,2},1),1)*[2 2 1 2], el{2,2}{1,2}]];
for i = 1 : 2
    if ~isempty(el{1,2}{i,i})
        consolidated_el = [consolidated_el; [ones(size(el{1,2}{i,i},1),1)*[1 2 i i], el{1,2}{i,i}]];
    end;
end;
consolidated_el = sortrows(consolidated_el, [7 -1 -2 -3 -4]);

while ~isempty(consolidated_el) & consolidated_el(1,7) < 0.2
    
    t1 = consolidated_el(1,1);
    t2 = consolidated_el(1,2);
    z1 = consolidated_el(1,3);
    z2 = consolidated_el(1,4);
    n1 = consolidated_el(1,5);
    n2 = consolidated_el(1,6);
    c{t1,t2}{z1,z2}(n1,n2) = 1;
    affected_edges = find_affected_edges(consolidated_el, 1, c);
    if ~isempty(affected_edges)
        affected_edges = intersect(consolidated_el(:,1:6), affected_edges, 'rows');
    end;
    c = update_corr_matrix(c, affected_edges);
    consolidated_el = remove_edges_from_consolidated_list(consolidated_el, affected_edges);
    rec_cost = consolidated_el(:,7);
    for i = 1:size(consolidated_el, 1)
        rec_cost(i) = recompute_labeling_cost(consolidated_el, i, c);
    end;
    consolidated_el(:,7) = rec_cost;
    consolidated_el = sortrows(consolidated_el, [7 -1 -3]);
end;

return;


