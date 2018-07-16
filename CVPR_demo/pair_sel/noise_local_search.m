function out = noise_local_search(lmat_struct2,np_array,num_cam)
%%%%% greedy_list_method
lf=[];
for j=1:length(lmat_struct2)
   lf=[lf;lmat_struct2{1,j}(:)]; 
    
end

ind=1:length(lf);
ind_greedy=zeros(size(lf));

pairs=nchoosek(1:length(np_array),2);
vertex_set=vertex_set_gen(np_array);

params_wadj.f=lf;
params_wadj.np_array=np_array;
params_wadj.pairs=pairs;
params_wadj.vertex_set=vertex_set;
params_wadj.ind=ind;


w_adjlist= weighted_adjacency_list_gen(params_wadj);

%% greedy iteration
for i=1:length(ind)

[u,v]=eval_end_vertices(ind(i),pairs,np_array,vertex_set);

u_adj=w_adjlist.vertex{1,u};
u_adj_w=w_adjlist.weight{1,u};

v_adj=w_adjlist.vertex{1,v};
v_adj_w=w_adjlist.weight{1,v};

common=intersect(u_adj,v_adj);
[row_u, ~]=find(u_adj'==common);
[row_v, ~]=find(v_adj'==common);

concat_w=[u_adj_w(row_u);v_adj_w(row_v)];
concat_w=concat_w(:,find(sum(concat_w>-1)));


dec_val=sum(concat_w(1,:).*concat_w(2,:));

if(dec_val>0)
    lf(i)=1;
elseif(dec_val<0)
    lf(i)=-1;
end
    

% if(i<3)
%  ind_greedy(i)=1;   
%  adjacency_list{u}=[adjacency_list{u} v];
%  adjacency_list{v}=[adjacency_list{v} u];
%    
% else
%     if(length(intersect(adjacency_list{u},adjacency_list{v}))<=0)
%       adjacency_list{u}=[adjacency_list{u} v];
%       adjacency_list{v}=[adjacency_list{v} u];
%       ind_greedy(i)=1;  
%     end
% 
% end



 
%%
 


    
    
    





end




out=pairwise_labelcollect(lf,np_array,num_cam);








end