function out =cell_normalizer_1d(cell1,k)
[a, b]=size(cell1);
out=cell(size(cell1));
for i=1:max(a,b)
out{i}=cell1{i}/k;
end







end