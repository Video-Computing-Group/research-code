function lcollect=pairwise_labelcollect(x_bip,np_array,num_cam)
pairs=nchoosek(1:num_cam,2);
bias=0;
for i=1:size(pairs,1)
 cam_a=np_array(pairs(i,1));
 cam_b=np_array(pairs(i,2));
 ct_np= cam_a*cam_b;
 lcollect{i}=sparse(reshape(x_bip(bias+1:bias+ct_np,1),cam_a,cam_b));
 bias=bias+ct_np;

    
end






end