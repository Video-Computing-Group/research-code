function out = max_prob_sel(f,lim)
%% max-prob selector 
[~,ind]=sort(f);
out=zeros(size(f));
out(ind(1:lim))=1;



end