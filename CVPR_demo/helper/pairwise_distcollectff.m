function [dist_collect, sim_collect]=pairwise_distcollectff(feat)
%%%%% [1-2,1-3,1-4,2-3........] 
camNamearr=1:size(feat,2);
pairs=nchoosek(camNamearr,2);
for i=1:size(pairs,1)
  cam_a=pairs(i,1);
  cam_b=pairs(i,2);
  dist_temp=sqdist(feat{cam_a}', feat{cam_b}',eye(size(feat{cam_a},2)));
  dist_collect{i}=dist_temp;
  sim_collect{i}=prob_mat(dist_temp,1*(median(dist_temp,2)));
  %mean(mean(dist_temp))
    
    
    
end





end