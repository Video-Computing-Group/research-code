function [pred_label,acc]= nn1clf_k(A,B,A_ID,B_ID,D)
pred_label=[];

if(nargin<5)
%[~,ind]=pdist2(B,A,'euclidean','smallest',1);
 D=eye(size(A,1));


end
Dmat=sqdist(A,B,D);
[~,inds]=sort(Dmat,2,'ascend');


temp_store=0;
rank_acc=zeros(size(inds));
for i=1:size(inds,1)
    temp_id_n=B_ID(inds(i,:));
    loc=min(find(temp_id_n==A_ID(i)));
    rank_acc(i,loc:end)=1;
    
end

acc=sum(rank_acc)/length(A_ID);

% for i=1:10
% 
% rank_acc(
% 
% pred_label=B_ID(ind);
% temp_ct=sum((A_ID==pred_label))+temp_store;
% 
% acc(i)=temp_ct/length(A_ID);
% temp_store=temp_ct;
% end



end
