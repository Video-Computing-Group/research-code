function sc = cform2SparseCell(pairwiseData,adjCell)

% convert the pairwise-data (like pairwise beliefs) written by c_inference
% to sparse_cell
N = size(adjCell,2);
sc = cell(N,N);
for i=1:N
    neighb_num = length(adjCell{i});
    for n=1:neighb_num
        j = adjCell{i}(n);
        if i<j
            data_ij = pairwiseData{i}{n};
            sc{i,j} = data_ij;
            sc{j,i} = data_ij';
        end
    end
end
