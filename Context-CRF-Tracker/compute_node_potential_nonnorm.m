function [local, states] = compute_node_potential_nonnorm(cord_1, cord_2, properties)

k = properties.Number_of_Nonzero_States;
if properties.Num_state_candidates > k
    K = properties.Num_state_candidates;
else 
    K = k;
end;
w = properties.NP_Shape_Distance_Weight;
lambda = properties.NP_Lambda;

local = cell(1, length(cord_1));
states = zeros(length(cord_1), k);
% q_seg = zeros(length(cord_1), 1);

% bp_map = double(imread(bp_im));
% bp_map = bp_map./255;

for c = 1:length(cord_1)
    np = zeros(k+1, 1);
    [cdts, f] = generate_candidates_feature_vector_prior(cord_1, c, cord_2, K);
    for cdt = 1:K
        if length(lambda) == 1
            np(cdt+1) = exp(-lambda*((1-w)*f(2*cdt-1)+w*f(2*cdt)));
        elseif length(lambda) > 1
            np(cdt+1) = exp((1-w)*-lambda(1)*f(2*cdt-1)+w*-lambda(2)*f(2*cdt));
        end;
    end; 
    l = [cdts', np(2:end)];
    l = sortrows(l, [-2]);
    np = [0, (l(1:k, 2))']';
    np(1) = max(0, 1-sum(np(2:end)));
%     q_seg(c) = compute_segmentation_quality(cord_1, c, bp_map, [1 0]);
%     np = q_seg(c) .* np;
    local{c} = np ./ sum(np);
    states(c, :) = (l(1:k, 1))';
end;

return;