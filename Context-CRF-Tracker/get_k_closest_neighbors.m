function [closest_nbor_cell_ids, distances] = get_k_closest_neighbors(cord_1, cell_1, cord_2, k)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% points_1 = cord_1(cell_1).PixelList;
% overlap = zeros(length(cord_2),1);
% distance = zeros(length(cord_2),1);
% 
% for i = 1:length(cord_2)
%     points_2 = cord_2(i).PixelList;
%     overlap(i) = size(intersect(fillInside(points_1), fillInside(points_2), 'rows'), 1);
%     distance(i) = norm(mean(fillInside(points_1), 1) - mean(fillInside(points_2), 1));
% end;
% 
% [o1, I1] = sort(overlap,1,'descend');
% I1 = I1(find(o1 > 0));
% closest_nbor_cell_ids = I1(1:min(3, length(I1)));
% if length(I1) < 3
%     distance(closest_nbor_cell_ids) = Inf;
%     [o2, I2] = sort(distance);
%     closest_nbor_cell_ids = [closest_nbor_cell_ids, I2(1: (3-length(I1)))];
% end;
% 
% return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

points_1 = cord_1(cell_1).PixelList;
distance = zeros(length(cord_2),3);

for i = 1:length(cord_2)
    points_2 = cord_2(i).PixelList;
    distance(i,1) = i;
    distance(i,2) = size(intersect(fillInside(points_1), fillInside(points_2), 'rows'), 1);
    distance(i,3) = norm(mean(fillInside(points_1), 1) - mean(fillInside(points_2), 1));
end;

distance = sortrows(distance, [-2, 3]);
closest_nbor_cell_ids = (distance(1:k, 1))';
distances = (distance(1:k, 2:3))';

return;