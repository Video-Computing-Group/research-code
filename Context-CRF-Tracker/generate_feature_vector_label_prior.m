function [feature, label] = generate_feature_vector_label_prior(cord_1, cell_1, cord_2, correspondence_1_2)

feature = [];
points_1 = cord_1(cell_1).PixelList;
h1 = compute_shape_histogram(shift_centroid(points_1, mean(fillInside(points_1), 1)), 'normalized_weighted_sector', [], [], -pi, pi, [], 8);

C = get_k_closest_neighbors(cord_1, cell_1, cord_2, 3);

for c = C'
    points_2 = cord_2(c).PixelList;                                                                                                                                                                                                                                                                                                              
    h2 = compute_shape_histogram(shift_centroid(points_2, mean(fillInside(points_2), 1)), 'normalized_weighted_sector', [], [], -pi, pi, [], 8);
    feature = [feature, norm(mean(fillInside(points_2), 1) - mean(fillInside(points_1), 1)), KLDiv(h2, h1)];
end;

if ~isempty(find(correspondence_1_2(cell_1, :) == 1))
    label = 1;
else
    label = -1;
end;

return;
 