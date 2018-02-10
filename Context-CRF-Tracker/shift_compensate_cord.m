function [cord_2, delta] = shift_compensate_cord(cord_1, cord_2, n)

points_1 = [];
points_2 = [];

figure(1), imshow(zeros(1024)); hold on;
for i = 1 : length(cord_1)
    plot(cord_1(i).PixelList(:,1), cord_1(i).PixelList(:,2), 'w.', 'MarkerSize', 2); hold on;
end;

figure(2), imshow(zeros(1024)); hold on;
for i = 1 : length(cord_2)
    plot(cord_2(i).PixelList(:,1), cord_2(i).PixelList(:,2), 'w.', 'MarkerSize', 2); hold on;
end;

for i = 1:n
    figure(1), points_1 = [points_1; ginput(1)];
    figure(2), points_2 = [points_2; ginput(1)];
end;

delta= mean(points_1 - points_2, 1);

for i = 1:length(cord_2)
    cord_2(i).PixelList = cord_2(i).PixelList + ones(length(cord_2(i).PixelList), 1)*delta;
end;

figure(3), imshow(zeros(1024)); hold on;
for i = 1 : length(cord_2)
    plot(cord_2(i).PixelList(:,1), cord_2(i).PixelList(:,2), 'w.', 'MarkerSize', 2); hold on;
end; 
return;



