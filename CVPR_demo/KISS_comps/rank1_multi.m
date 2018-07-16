function [ accuracy_r1 ] = rank1_multi(M, data,idxa, idxb,ind_a,ind_b)


M=eye(size(data,1));

dist=Distmat((data(:,idxa))', (data(:,idxb))',M);
[~,inds]=sort(dist,2,'ascend');
inds=inds(:,1);
result=0;

for i=1:size(dist,1)
if (ind_a((i))==ind_b(inds(i)))
    
  result=result+1
    
end
    
    
end
accuracy_r1=100*(result/i);






end


