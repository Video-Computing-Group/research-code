function [boundary] = get2DBoundary(points, method)

% The method can be either of 'union' or 'intersection'
x = points(:,1);
y = points(:,2);

x_u = unique(points(:,1));
y_u = unique(points(:,2));

boundary_x = [];
for i = 1:length(x_u)
    yy = [min(y(find(x == x_u(i)))) max(y(find(x == x_u(i))))]';
    xx = x_u(i).*ones(length(yy),1);
    boundary_x = [boundary_x; [xx, yy]];
end;

boundary_y = [];
for i = 1:length(y_u)
    xx = [min(x(find(y == y_u(i)))) max(x(find(y == y_u(i))))]';
    yy = y_u(i).*ones(length(xx),1);
    boundary_y = [boundary_y; [xx, yy]];
end;

if strcmpi(method, 'intersection')
    boundary = intersect(boundary_x, boundary_y, 'rows');
elseif strcmpi(method, 'union')
    boundary = union(boundary_x, boundary_y, 'rows'); 
else
    disp('Method could be either union or intersection');
end;

return;