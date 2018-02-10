function [np, states] = compute_node_potential_ptive_only(cord_1, cord_2, feature_set_positive, k)

np = zeros(length(cord_1), k);
states = zeros(length(cord_1), k);

for c = 1:length(cord_1)
    [cdts, f] = generate_candidates_feature_vector_prior(cord_1, c, cord_2, k);
    states(c, :) = cdts;
    p_p = zeros(1,k);
    for cdt = 1:k
        %p_p(cdt) = 1/mahal(f(2*cdt-1:2*cdt), feature_set_positive(:,3:4));
        p_p(cdt) = exp(-f(2*cdt-1)) ;%+ 0.7*exp(-f(2*cdt));
    end;
    for cdt = 1:k
        np(c, cdt) = p_p(cdt)/sum(p_p);
    end;
end;

return;