function [local, states] = compute_node_potential_nonnorm_temporal(cord_1, cord_2, properties)

k = properties.Number_of_Nonzero_States;
w = properties.NP_Shape_Distance_Weight;
lambda = properties.NP_Lambda;

local = cell(1, length(cord_1));
states = zeros(length(cord_1), k);
% q_seg = zeros(length(cord_1), 1);

% bp_map = double(imread(bp_im));
% bp_map = bp_map./255;

for c = 1:length(cord_1)
    np = zeros(k, 1);
    [cdts, f] = generate_candidates_feature_vector_prior(cord_1, c, cord_2, k);
    states(c, :) = cdts;
    for cdt = 1:k
        if length(lambda) == 1
            np(cdt) = exp(-lambda*((1-w)*f(2*cdt-1)+w*f(2*cdt)));
        elseif length(lambda) > 1
            np(cdt) = exp((1-w)*-lambda(1)*f(2*cdt-1)+w*-lambda(2)*f(2*cdt));
        end;
    end;
    %     q_seg(c) = compute_segmentation_quality(cord_1, c, bp_map, [1 0]);
    %     np = q_seg(c) .* np;
    local{c} = np;
end;

return;