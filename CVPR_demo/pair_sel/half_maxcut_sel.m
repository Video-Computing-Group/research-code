function ind_half = half_maxcut_sel(f,np_array,lim)
%%% 1/2 approximation algo. for triangle free subset selection

[f_sorted,ind]=sort(f);
f_sorted(lim+1:end,:)=0;
ind_half=zeros(size(f));

pairs=nchoosek(1:length(np_array),2);
vertex_set=vertex_set_gen(np_array);

params_wadj.f=f_sorted;
params_wadj.np_array=np_array;
params_wadj.pairs=pairs;
params_wadj.vertex_set=vertex_set;
params_wadj.ind=ind;


w_adjlist= weighted_adjacency_list_gen(params_wadj);


set_a=1;
set_b=2;
for i=3:sum(np_array)
if (size(w_adjlist.vertex{i},1)>0)
 [~,wt_indices_b]=find(w_adjlist.vertex{i}==set_a');   
 [~,wt_indices_a]=find(w_adjlist.vertex{i}==set_b');
 
 wt_b=w_adjlist.weight{i}(wt_indices_b);
 wt_a=w_adjlist.weight{i}(wt_indices_a);
 
 if(sum(wt_b)<sum(wt_a))
    set_b=[set_b i];
    ind_half(w_adjlist.edge{i}(wt_indices_b))=1;
 else
    set_a=[set_a i]; 
    ind_half(w_adjlist.edge{i}(wt_indices_a))=1;

 end
end

    
    
end





end