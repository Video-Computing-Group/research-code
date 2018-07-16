function [tot_q,chkmt]=QAT(D,L1,Indq)



tot_q=0;
[~,Ind]=sort(D,2,'ascend');
 for i=1:length(Indq)
 i;
 ct_q=find(Indq(Ind(i,:))==L1((i)));
 tot_q=tot_q+ct_q;   
 chkmt(1,i)=Indq(i);
 chkmt(2,i)=ct_q;
 end

end