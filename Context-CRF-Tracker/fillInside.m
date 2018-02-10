function filled_Points = fillInside(points)

x = points(:,1);
y = points(:,2);

x_u = unique([min(points(:,1)):0.1:max(points(:,1))]);
y_u = unique([min(points(:,2)):0.1:max(points(:,2))]);

% x_u = unique(points(:,1));
% y_u = unique(points(:,2));

filled_Points_x = [];

for i = 1:length(x_u)
    
    yy = [min(y(find(x == x_u(i)))):max(y(find(x == x_u(i))))]';
    xx = x_u(i).*ones(length(yy),1);
    
    filled_Points_x = [filled_Points_x; [xx, yy]];
    
end;


filled_Points_y = [];

for i = 1:length(y_u)
    
    xx = [min(x(find(y == y_u(i)))):max(x(find(y == y_u(i))))]';
    yy = y_u(i).*ones(length(xx),1);
    
    filled_Points_y = [filled_Points_y; [xx, yy]];
    
end;

filled_Points = intersect(filled_Points_x, filled_Points_y, 'rows');

return;