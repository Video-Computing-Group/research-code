function [adj] = adjacency_4m_timestamps(obslist)

adj = zeros(size(obslist, 1), size(obslist, 1));
for i = 1 : size(obslist, 1)
    for j = 1 : size(obslist, 1)
        ti_range = [obslist(i,2):obslist(i,3)];
        tj_range = [obslist(j,2):obslist(j,3)];
        t_comm = intersect(ti_range, tj_range);
        if ~isempty(t_comm) && i~=j %&& (length(t_comm)/length(ti_range)) >=0.2 && (length(t_comm)/length(tj_range)) >=0.2
            adj(i,j) = 1;
        end;
    end;
end;

return;