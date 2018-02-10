function [prior_dist] = compute_prior_dist(cord_1, cord_2, cell_1)

prior_dist = zeros(1, length(cord_2));
points_1 = cord_1(cell_1).PixelList;
weight = 1;

for i = 1:length(cord_2)
    points_2 = cord_2(i).PixelList;
    overlap = size(intersect(fillInside(points_1), fillInside(points_2), 'rows'), 1);
    distance_location = norm(mean(fillInside(points_1), 1) - mean(fillInside(points_2), 1));
    overlap = min(overlap/size(fillInside(points_1),1), overlap/size(fillInside(points_2),1));
    if overlap > 0 | distance_location < 35
         %distance_shape = ModHausdorffDist(points_1, points_2);
         distance_shape = ModHausdorffDist(shift_centroid(points_1, mean(fillInside(points_1), 1)), ...
                                           shift_centroid(points_2, mean(fillInside(points_2), 1)));
         r = cord_1(cell_1).EquivDiameter / 2;
         prior_dist(i) = weight*exp(-distance_shape/r) + (1-weight)*overlap;
%         distance_shape = ModHausdorffDist(shift_centroid(points_1, mean(fillInside(points_1), 1)), ...
%                                           shift_centroid(points_2, mean(fillInside(points_2), 1)));
%         prior_dist(i) = weight*exp(-distance_shape) + (1-weight)*exp(-distance_location);
%         prior_dist(i) = exp(-distance_shape)*exp(-distance_location);
%         ov = min(overlap/size(fillInside(points_1),1), overlap/size(fillInside(points_2),1));
%         prior_dist(i) = exp((1-ov)*distance_shape);
%         [R,T] = icp(points_1, points_2);
%         prior_dist(i) = exp(-norm(T)) + exp(-norm(R-eye(2), 'fro'));
    end;
end;

prior_dist = prior_dist ./ sum(prior_dist);
return;