function [candidates, features] = generate_candidates_feature_vector_prior(cord_1, cell_1, cord_2, k)

features = [];
points_1 = cord_1(cell_1).PixelList;
h1 = compute_shape_histogram(shift_centroid(points_1, mean(fillInside(points_1), 1)), 'normalized_weighted_sector', [], [], -pi, pi, [], 8);

C = get_k_closest_neighbors(cord_1, cell_1, cord_2, k);
candidates = C;

for c = C
    points_2 = cord_2(c).PixelList;                                                                                                                                                                                                                                                                                                              
    h2 = compute_shape_histogram(shift_centroid(points_2, mean(fillInside(points_2), 1)), 'normalized_weighted_sector', [], [], -pi, pi, [], 8);
    features = [features, norm(mean(fillInside(points_2), 1) - mean(fillInside(points_1), 1)), KLDiv(h2, h1)];
end;

return;