function [pred_label,acc]= nn1clf(A,B,A_ID,B_ID,D)

if(nargin<5)
[~,ind]=pdist2(B,A,'euclidean','smallest',1);
else
Dmat=Distmat(A,B,D);
[~,inds]=sort(Dmat,2,'ascend');
ind=inds(:,1);
end
pred_label=B_ID(ind);
acc=sum((A_ID==pred_label))/length(A_ID);





end