function [isdivided, d1, d2] = detect_cell_divisions(cord_1, parent, cord_2, children_pair, threshold)

isdivided = 0;

c1 = parent;
c11 = children_pair(1);
c12 = children_pair(2);

%% MHD computation between boundaries of the paernt and the combined
%% boundaries of the children pair.

comb_pts = union(cord_2(c11).PixelList, cord_2(c12).PixelList, 'rows');
comb_boundary = get2DBoundary(comb_pts, 'union');
d1 = ModHausdorffDist(shift_centroid(comb_boundary, mean(fillInside(comb_boundary), 1)), shift_centroid(cord_1(c1).PixelList, mean(fillInside(cord_1(c1).PixelList), 1)));


%% Distance computed using area-ratios.

d2 = 0.5*(abs(0.5 - cord_2(c11).Area/cord_1(c1).Area) + abs(0.5 - cord_2(c12).Area/cord_1(c1).Area));

%% Check if the combined distance is less than a threshold.

if d1 <= threshold(1) & d2 <= threshold(2)
    isdivided = 1;
end;

return;