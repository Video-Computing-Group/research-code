function [local, states] = compute_node_potential(cord_1, cord_2, feature_set_positive, feature_set_negative, k)

local = cell(1, length(cord_1));
states = zeros(length(cord_1), k);
% q_seg = zeros(length(cord_1), 1);

% bp_map = double(imread(bp_im));
% bp_map = bp_map./255;

for c = 1:length(cord_1)
    np = zeros(k+1, 1);
    [cdts, f] = generate_candidates_feature_vector_prior(cord_1, c, cord_2, k);
    states(c, :) = cdts;
    p_p = zeros(1,k);
    p_n = zeros(1,k);
    for cdt = 1:k
        d_p = mahal(f(2*cdt-1:2*cdt), feature_set_positive(:,3:4));
        d_n = mahal(f(2*cdt-1:2*cdt), feature_set_negative(:,3:4));
        p_p(cdt) = d_n/(d_p+d_n);
        p_n(cdt) = d_p/(d_p+d_n);
    end;
    np(1) = min(p_n);
    np_rest = 1 - np(1);
    for cdt = 1:k
        np(cdt+1) = p_p(cdt)*np_rest/sum(p_p);
    end;
%     q_seg(c) = compute_segmentation_quality(cord_1, c, bp_map, [1 0]);
%     np = q_seg(c) .* np;
    local{c} = np;
end;

return;
        
    