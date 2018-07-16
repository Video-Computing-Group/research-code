function out = overlap_stat_gen(cam_ind)
pairs=nchoosek(1:length(cam_ind),2);
triplets=nchoosek(1:length(cam_ind),3);
pairwise_overlap=0;
triplewise_overlap=0;

for i=1:size(pairs,1)
 temp1=intersect(cam_ind{pairs(i,1)},cam_ind{pairs(i,2)});
 pairwise_overlap=pairwise_overlap+length(temp1)/max(length(cam_ind{pairs(i,1)}),length(cam_ind{pairs(i,2)}));
    
end
pairwise_overlap_mean=100*(pairwise_overlap/size(pairs,1));

for j=1:size(triplets,1)
 temp1=intersect(cam_ind{triplets(j,1)},cam_ind{triplets(j,2)});
 temp2=intersect(cam_ind{triplets(j,3)},temp1);
 triplewise_overlap=triplewise_overlap+length(temp2)/max([length(cam_ind{triplets(j,1)}),length(cam_ind{triplets(j,2)}),length(cam_ind{triplets(j,3)})]);

 
end
triplewise_overlap_mean=100*(triplewise_overlap/size(triplets,1));

out.pair_stat=pairwise_overlap_mean;
out.triple_stat=triplewise_overlap_mean;

end