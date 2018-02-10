function [feature_set] = compute_features_prior(points_1, points_2)

 feature_set = zeros(1, 4);
 feature_set(1) = size(intersect(fillInside(points_1), fillInside(points_2), 'rows'), 1) ...
                  / size(fillInside(points_1), 1);
 feature_set(2) = size(intersect(fillInside(points_1), fillInside(points_2), 'rows'), 1) ...
                  / size(fillInside(points_2), 1);
 feature_set(3) = norm(mean(fillInside(points_1), 1) - mean(fillInside(points_2), 1));
 feature_set(4) = ModHausdorffDist(shift_centroid(points_1, mean(fillInside(points_1), 1)), shift_centroid(points_2, mean(fillInside(points_2), 1)));

%feature_set = ModHausdorffDist(points_1, points_2);

%feature_set(4) = compute_l2_distance_GMM_shape(points_1, points_2, 3);

return;