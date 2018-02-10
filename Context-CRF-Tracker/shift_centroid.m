function [points_shifted] = shift_centroid(points, new_centroid)

points_shifted = points - ones(size(points,1),1)*new_centroid;
return;