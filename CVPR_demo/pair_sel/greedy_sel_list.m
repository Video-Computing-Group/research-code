function out = greedy_sel_list(f,np_array,lim)
%%%%% greedy_list_method


[~,ind]=sort(f);
ind_greedy=zeros(size(f));

pairs=nchoosek(1:length(np_array),2);
vertex_set=vertex_set_gen(np_array);
adjacency_list=cell(1,sum(np_array));

%% greedy iteration
i=1;
while sum(ind_greedy)<lim

[u,v]=eval_end_vertices(ind(i),pairs,np_array,vertex_set);

if(i<3)
 ind_greedy(i)=1;   
 adjacency_list{u}=[adjacency_list{u} v];
 adjacency_list{v}=[adjacency_list{v} u];
   
else
    if(length(intersect(adjacency_list{u},adjacency_list{v}))<=0)
      adjacency_list{u}=[adjacency_list{u} v];
      adjacency_list{v}=[adjacency_list{v} u];
      ind_greedy(i)=1;  
    end

end



 
%%
 i=i+1; %%incrementing counter
 
 %% breaking while loop
 if(i>length(f))
     break
 end
 %%
 
    
    
    





end




out=zeros(size(f));
out(ind(find(ind_greedy)))=1;









end