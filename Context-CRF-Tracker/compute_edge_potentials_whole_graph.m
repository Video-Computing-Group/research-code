function [ep] = compute_edge_potentials_whole_graph(cord_1, cord_2, states, adj, properties, mode)

gamma = properties.EP_Lambda;

ep = cell(length(cord_1), length(cord_1));

for i = 1:size(ep, 1)
    J = find(adj(i,:) == 1);
    for j = J
        if adj(i, j) == 1
            states_1 = states(i, :);
            states_2 = states(j, :);
            if strcmpi(mode, 'spatial')
                ep{i, j} = compute_edge_potential_gamma(cord_1, cord_2, i, j, states_1, states_2, gamma);
            elseif strcmpi(mode, 'temporal')
                ep{i, j} = compute_edge_potential_gamma_temporal(cord_1, cord_2, i, j, states_1, states_2, gamma);
            end;
        end;
    end;
end;

return;