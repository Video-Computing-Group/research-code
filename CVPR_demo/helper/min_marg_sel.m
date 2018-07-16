function out = min_marg_sel(sim_struct,num_pairs,lim)
%% min margin select
f=[];
for i=1:num_pairs
    temp_f=abs(sim_struct{i}-(1-sim_struct{i}));
    f=[f;temp_f(:) ];
    
end

[~,ind]=sort(f);
out=zeros(size(f));
out(ind(1:lim))=1;




end