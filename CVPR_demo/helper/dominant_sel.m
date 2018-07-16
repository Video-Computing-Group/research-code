

function out = dominant_sel(Feat,num_cam,lim)
%%%% dominant clustering based selector


pairs=nchoosek(1:num_cam,2);
lim_per_pair=round(lim/size(pairs,1));
out=[];
for i=1:size(pairs,1)
i
per_probe_lim=lim_per_pair/size(Feat{pairs(i,1)},1);   
pair_wise_mat_sel=zeros(size(Feat{pairs(i,1)},1),    size(Feat{pairs(i,2)},1));
 for j=1:size(Feat{pairs(i,1)},1)
     
     Temp_mat=[Feat{pairs(i,1)}(j,:);Feat{pairs(i,2)}];
     
     [S] = create_sim_matrix(Temp_mat',1);
     [clusters ,~ ,~,~,ncluster] = clusterDS(S, 'MaxClust', 2);
     clusters_temp=clusters(2:end);
     [ind_temp]=find(clusters_temp==1);
     pair_wise_mat_sel(j,ind_temp(1:per_probe_lim))=1;
     
 end   
    
 out=[out; pair_wise_mat_sel(:)];
    
    

    
    
    
    
    
    
    
    
    
    
    
    
end



end