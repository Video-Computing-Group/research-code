function lmat =lmat_gen(cam_a_id,cam_b_id)
%label matrix generator 
lmat=sparse(length(cam_a_id),length(cam_b_id)); %% was zeros before
for i=1:length(cam_a_id)
    ind=find(cam_b_id==cam_a_id(i));
    lmat(i,ind)=1;








end