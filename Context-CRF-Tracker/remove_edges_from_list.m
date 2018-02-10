function [el] = remove_edges_from_list(el, n1, n2)

i = find(el(:,1)==n1 | el(:,2)==n2);
el(i,:) = [];
el = sortrows(el, [3 -1]);

return;