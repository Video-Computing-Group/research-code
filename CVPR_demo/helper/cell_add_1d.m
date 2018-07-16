function out =cell_add_1d(cell1,cell2)
[a, b]=size(cell1);
out=cell(size(cell1));
for i=1:max(a,b)
if(size(cell1{i},1)>0 && size(cell1{i},2)>0)
out{i}=cell1{i}+cell2{i};
else
out{i}=zeros(size(cell2{i}))+cell2{i};    
end


end




end