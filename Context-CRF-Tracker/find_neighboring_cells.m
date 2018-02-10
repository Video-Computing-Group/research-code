function [neighbors] = find_neighboring_cells(cord, cell_id, method)

warning off;
neighbors = [];

if strcmpi(method, 'ellipse')
    [A, c] = MinVolEllipse((cord(cell_id).PixelList)', 0.0001);
elseif strcmpi(method, 'circle')
    [c, R]  = minboundcircle(cord(cell_id).PixelList(:,1), cord(cell_id).PixelList(:,2));
    c = c';
    A = (1/R^2).*eye(2);
else
    disp('Please enter a valid method name: ''ellipse'' or ''circle''');
end;

for i = 1:length(cord)
    D = mhd_computation_pt_cloud_single_cell(fillInside(cord(i).PixelList), A, c);
    if length(find(D<=1)) > 0
        neighbors = [neighbors, i];
    end;
end;

neighbors(find(neighbors==cell_id)) = [];

return;