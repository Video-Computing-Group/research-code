function [u,v]=eval_end_vertices(k,pairs,np_array,vertex_set)
%%% gives vertex indices from k-th edge
pwise_edgect=edge_count_pairwise(pairs,np_array);
pwise_cmedgect=cumsum(pwise_edgect);
for i=1:length(pwise_cmedgect)
    if(k<=pwise_cmedgect(i))
        cam_pair_id=i;
        if(i>1)
        pair_specific_lin_cam_pair_ind=k-pwise_cmedgect(i-1);
        else
        pair_specific_lin_cam_pair_ind=k;
        end
        break;
    end
    
end
cam_a=pairs(cam_pair_id,1);
cam_b=pairs(cam_pair_id,2);
sz1=np_array(cam_a);
sz2=np_array(cam_b);


[u_temp,v_temp]=ind2sub([sz1 sz2],pair_specific_lin_cam_pair_ind);
u=vertex_set{cam_a}(u_temp);
v=vertex_set{cam_b}(v_temp);
end