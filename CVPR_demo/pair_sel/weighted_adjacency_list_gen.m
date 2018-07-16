function wadjlist= weighted_adjacency_list_gen(params)
%% creates weighted adjacency list from edge weight list f and vertex set 
%%



np_array=params.np_array;
pairs=params.pairs;
vertex_set=params.vertex_set;
f=params.f;
ind=params.ind;



 wadjlist.vertex=cell(1,sum(np_array));
 wadjlist.weight=cell(1,sum(np_array));
 wadjlist.edge=cell(1,sum(np_array));

for i=1:length(f)
    
% if(ind(i)==1435)
%     pause
% end
if(abs(f(i))>0) %%was >0 before noise
[u,v]=eval_end_vertices(ind(i),pairs,np_array,vertex_set);


 wadjlist.vertex{u}=[ wadjlist.vertex{u} v];
 wadjlist.vertex{v}=[ wadjlist.vertex{v} u];
 
    
 wadjlist.weight{u}=[wadjlist.weight{u} f(i)];
 wadjlist.weight{v}=[wadjlist.weight{v} f(i)];
 
 
 wadjlist.edge{u}=[ wadjlist.edge{u} ind(i)];
 wadjlist.edge{v}=[ wadjlist.edge{v} ind(i)];
end
end



end