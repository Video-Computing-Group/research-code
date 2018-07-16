%%% performs a greedy max-cut with at least 1/2 approximation guarantee and
%%% select edges for labelling 
%% 

function [label_temp,cam_vind_seta,cam_vind_setb] = greedy_maxcut_sel(w,np_array)
% rng('default')
pairs=nchoosek(1:length(np_array),2);
pairs_merged=str2num(strcat(num2str(pairs(:,1)),num2str(pairs(:,2))));
per_pair_2ndcam=cumsum(np_array(pairs(:,2)));

num_v=sum(np_array);
vertex_ind=1:num_v;
cam_vind=repelem(1:length(np_array),np_array); %%% useful function
camwise_vpind=vertex_ind-repelem(cumsum([0 np_array(2:end)]),np_array); %% corresponding loop-version of this line also runs with same speed 
label_temp=zeros(size(w));
a_cnt=1;
b_cnt=1;
for i=1:num_v
random_vsel=randi(num_v-i+1);
random_vind=vertex_ind(random_vsel);
vertex_ind(random_vsel)=[];

if(i==1)
    set_a(a_cnt)=random_vind;
    a_cnt=a_cnt+1;
end
if(i==2)
    set_b(b_cnt)=random_vind;
    b_cnt=b_cnt+1;

end

if(i>=3)
    
%%% ct cameras
temp_cama=cam_vind(set_a);
temp_camb=cam_vind(set_b);
temp_camcta=repmat(cam_vind(random_vind),size(temp_cama));
temp_camctb=repmat(cam_vind(random_vind),size(temp_camb));

%%% persons


temp_camap=camwise_vpind(set_a);
temp_cambp=camwise_vpind(set_b);
temp_camctpa=repmat(camwise_vpind(random_vind),size(temp_camap));
temp_camctpb=repmat(camwise_vpind(random_vind),size(temp_cambp));

%%% same-partite deletion 
partite_samea=find(temp_cama==temp_camcta);
partite_sameb=find(temp_camb==temp_camctb);
temp_cama(partite_samea)=[];%cam deletion 
temp_camb(partite_sameb)=[];%cam deletion

temp_camcta(partite_samea)=[];%cam deletion 
temp_camctb(partite_sameb)=[];%cam deletion


temp_camctpa(partite_samea)=[];%ct person deletion 
temp_camctpb(partite_sameb)=[];%ct person deletion

temp_camap(partite_samea)=[];%%person deletion
temp_cambp(partite_sameb)=[];%%person deletion 
% sm_cama(partite_samea)=[];
% sm_camb(partite_sameb)=[];


 dsm_cama=find(temp_cama>temp_camcta);
 dsm_camb=find(temp_camb>temp_camctb);
 
if(~isempty(temp_cama))
cam_paira=str2num(strcat(num2str(min(temp_cama',temp_camcta(1))),num2str(max(temp_cama',temp_camcta(1)))));

[shift_a,~]=find(pairs_merged==cam_paira');
[temp_camap(dsm_cama) temp_camctpa(dsm_cama)] = deal(temp_camctpa(dsm_cama),temp_camap(dsm_cama) );%% swap elements
temp_camctpa=temp_camctpa+np_array(1)*(shift_a'-1);
edge_inda=sub2ind([np_array(1) sum(np_array(pairs(:,2))) ], temp_camap', temp_camctpa');
else
edge_inda=[];
end

if(~isempty(temp_camb))
cam_pairb=str2num(strcat(num2str(min(temp_camb',temp_camctb(1))),num2str(max(temp_camb',temp_camctb(1)))));
[shift_b,~]=find(pairs_merged==cam_pairb');
[temp_cambp(dsm_camb) temp_camctpb(dsm_camb)] = deal(temp_camctpb(dsm_camb),temp_cambp(dsm_camb) );%% swap elements
% i
temp_camctpb=temp_camctpb+np_array(1)*(shift_b'-1);
edge_indb=sub2ind([np_array(1) sum(np_array(pairs(:,2)))  ], temp_cambp', temp_camctpb');
else
edge_indb=[];    
end






if(sum(w(edge_inda))>sum(w(edge_indb)))
 label_temp(edge_inda)=1;
 set_b(b_cnt)=random_vind;
 b_cnt=b_cnt+1;
%  wt=wt+sum(w(edge_inda))
%  w'*label_temp
else
 label_temp(edge_indb)=1;
 set_a(a_cnt)=random_vind;
 a_cnt=a_cnt+1;

 

    
end    
    
end

end


%%

%%



cam_vind_seta=cam_vind(set_a);
cam_vind_setb=cam_vind(set_b);






    
end





    
    
    







