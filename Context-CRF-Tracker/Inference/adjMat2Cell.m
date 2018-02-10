function adjCell = adjMat2Cell(adjMat)

% create the cell-array of neighbours indices
N = size(adjMat,1);
adjCell = cell(1,N);
for i=1:N
    adjCell{i} = find(adjMat(i,:));
end
