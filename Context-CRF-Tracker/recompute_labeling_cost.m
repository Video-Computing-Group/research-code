function [rec_cost] = recompute_labeling_cost(consolidated_el, ednum, correspondence_mat)

affected_edges = find_affected_edges(consolidated_el, ednum, correspondence_mat);

if isempty(affected_edges)
    rec_cost = consolidated_el(ednum, 7);
else
    [~, rownum, ~] = intersect(consolidated_el(:,1:6), affected_edges, 'rows');
    extra_costs = consolidated_el(rownum, 7);
    rec_cost = (consolidated_el(ednum, 7) + sum(extra_costs))/(length(extra_costs)+1);
end;

return;