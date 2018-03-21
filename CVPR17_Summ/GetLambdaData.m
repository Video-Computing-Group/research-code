%--------------------------------------------------------------------------
% This code takes a DxN matrix of N data points in a D-dimensional 
% space and returns the regularization constant of the data_21 norm
% Y: DxN data matrix
% lambda1: regularization parameter for data_21 norm
%-------------------------------------------------------------------------
%--------------------------------------------------------------------------

function lambda = GetLambdaData(Y)

[~,N] = size(Y);
T = zeros(N,1);
for i = 1:N
    yi = Y(:,i);
    T(i) = norm(yi' * Y);
end
lambda = max(T);

