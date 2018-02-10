function [points_shifted] = shift_centroid_to_mean(points)

new_centroid = mean(fillInside(points), 1);
points_shifted = shift_centroid(points, new_centroid);
return;