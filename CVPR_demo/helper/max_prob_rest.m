function out2 = max_prob_rest(f,x_bip,lim1,lim2)
%% max-prob selector 
if(size(f,1)~=size(x_bip,1))
    x_bip=x_bip';
end
f=f.*not(x_bip);

[~,ind]=sort(f);
out=zeros(size(f));
out(ind(1:lim2-lim1))=1;

out2=x_bip+out;




end