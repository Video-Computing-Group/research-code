function [Idx_A ,Idx_B,La]= completeind2ff(asmat,lmat,idxtrain1,idxtrain2,idxa,idxb)
count0=1;
count1=1;
% Idx_A1=[];
% Idx_A0=[];
% Idx_B1=[];
% Idx_B0=[];
% La1=[];
% La0=[];
% for i =1:size(asmat,1)
%     for j=1:size(lmat,2)
%         if(asmat(i,j)==1 && lmat(i,j)==0)
%         Idx_A0(1,count0)=idxa(idxtrain1(i));

asmat=full(asmat);
lmat=full(lmat);
new2=asmat(:).*lmat(:);
[posind]=find(new2==1);
[posi, posj]=ind2sub(size(lmat),posind');
Idx_A1=idxa(idxtrain1(posi));
Idx_B1=idxb(idxtrain2(posj));
La1=ones(1,length(posi));

if(sum(sum(isnan(lmat)))>0)
temp_xx=find(lmat(:)==1);
temp_var=zeros(size(lmat(:)));
temp_var(temp_xx)=1;
temp_var=not(temp_var);

new1=asmat(:).*temp_var;
else
new1=asmat(:).*not(lmat(:));    
    
end
[negind]=find(new1==1);
[negi, negj]=ind2sub(size(lmat),negind');
Idx_A0=idxa(idxtrain1(negi));
Idx_B0=idxb(idxtrain2(negj));
La0=zeros(1,length(negi));




 
Idx_A=[Idx_A1 Idx_A0];
Idx_B=[Idx_B1 Idx_B0];
La=[La1 La0];
La=logical(La);








end