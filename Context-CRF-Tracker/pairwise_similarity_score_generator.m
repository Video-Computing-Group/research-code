function [corr_mat, edge_list, div_list] = pairwise_similarity_score_generator(cord_1, cord_2, mode, properties)

corr_mat = zeros(length(cord_1), length(cord_2));
cells_1 = [1:length(cord_1)]';
cells_2 = [1:length(cord_2)]';
if strcmpi(mode, 'temporal')
    disp('Starting Division Event Detection ...');
    [division_results] = enlist_cell_division_events(cord_1, cord_2, properties);
    if ~isempty(division_results)
        a1 = division_results(:,1);
        a2 = division_results(:,2:3);
        cells_1(a1(:)) = [];
        cord_1(a1(:)) = [];
        cells_2(a2(:)) = [];
        cord_2(a2(:)) = [];
    end;
    disp('Division Event Detection Complete ...');
    div_list = [];   % New
    for i = 1:size(division_results,1)
        n1 = division_results(i,1);
        for n2 = division_results(i,2:3)
            %corr_mat(n1, n2) = 1;
            corr_mat(n1, n2) = 0; % Changed last line
            div_list = [div_list; [n1 n2]]; % New
        end;
    end;
end;

[adj] = compute_adjacency_matrix(cord_1);

disp('Computing Node and Edge Potentials ...');
if strcmpi(mode, 'spatial')
    [np, states] = compute_node_potential_nonnorm(cord_1, cord_2, properties);
elseif strcmpi(mode, 'temporal')
    [np, states] = compute_node_potential_nonnorm(cord_1, cord_2, properties);
end;

[ep] = compute_edge_potentials_whole_graph(cord_1, cord_2, states, adj, properties, 'spatial');
disp('Node and Edge Potentials Computation Complete ...');

disp('Starting Inference ...');
[bel, converged, belpairs, msg] = inference(adj, ep, np, 'loopy', 'temperature', properties.MRF_Temperature, ...
    'sum_or_max', properties.MRF_Sum_or_Max, 'strategy', properties.MRF_Strategy, ...
    'log', properties.MRF_Log, 'log_bels', properties.MRF_Log_Bels);
disp('Inference Complete ...');

edge_list = zeros(size(states,1)*size(states,2), 2);
counter = 1;
for n1 = 1 : size(states,1)
    for n2 = states(n1, :)
        edge_list(counter, 1) = cells_1(n1);
        edge_list(counter, 2) = cells_2(n2);
        edge_list(counter, 3) = bel{n1}(find(states(n1,:)==n2)+1);
        counter = counter + 1;
    end;
end;
edge_list = sortrows(edge_list, [-3 1]);

return;
