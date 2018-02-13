% Andreas Krause (krausea@gmail.com)
% Compute Euclidean distances between a set of vectors
%
% function D = sfo_dist(coords)
% coords: N x D matrix of points (each of the N rows is a vector)
% returns an N x N matrix of Euclidean distances
%
% Example: coords = rand(10,2); D = sfo_dist(coords)

function D = sfo_dist(coords)
n = size(coords,1);
D = zeros(n);

for i = 1:n
    x = coords(i,:);
    D(i,:) = sqrt(sum((ones(n,1)*x-coords).^2,2));
end
