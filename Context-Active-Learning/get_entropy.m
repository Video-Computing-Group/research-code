function H = get_entropy(X)

% Add a jitter for numerical stability
% X = X + 1e-6;

% Establish size of data
m = size(X,2);

% Housekeeping
H = zeros(1,m);

for Column = 1:m,
    P = X(:,Column);
    % P = sort(P,'descend');
    % P = P(1:3);
    temp = log2(P);
    temp(isinf(temp)) = 0;
    temp(isnan(temp)) = 0;
    H(Column) = -sum(P .* temp);
end

H(isnan(H)) = 0;

