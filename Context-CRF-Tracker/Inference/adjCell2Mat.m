function adjMat = adjCell2Mat(adjCell)

% create a sparse adjacencies matrix from the cell-array of neighbours indices
N = length(adjCell);
adjMat = sparse(N,N);
for i=1:N
    adjMat(i,adjCell{i}) = 1;
end
