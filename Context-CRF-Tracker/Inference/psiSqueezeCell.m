function psiCell = psiSqueezeCell(lambda,adjCell,potts_model)

% squeeze the psi-data into the form read by c_inference
if (nargin<3)
    potts_model = 0;
end
N = size(adjCell,2);
psiCell = cell(1,N);
for i=1:N
    if potts_model
        psiCell{i} = lambda(i,adjCell{i});
    else        
        neighb_num = length(adjCell{i});
        psiCell{i} = cell(1,neighb_num);
        for n=1:neighb_num
            j = adjCell{i}(n);
            psiCell{i}{n} = lambda{i,j};
        end
    end
end
