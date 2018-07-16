%%%%%% 
function [ACC_mat,Req_vec,label_gain_vec,stats] = Ward_sim(params)
load Ward_dataset_folds3.mat
Train_mats=Ward_dataset_folds2.Train_single{1,params.num_folds};
Train_matm=Ward_dataset_folds2.Train_multi{1,params.num_folds};
%Test_mats=Raid_dataset_folds.Test_single{1,params.num_folds};
Test_matm=Ward_dataset_folds2.Test_multi{1,params.num_folds};

% Test_matm=Test_mats;
%%
% if params.randomize_split==1
%  total_id=[Train_mats.camp{1} Test_mats.camp{1}];
%  random_index=randperm(length(total_id));
%  total_id=total_id(random_index);
%  train_rand_id=total_id(1:25);
%  test_rand_id=total_id(26:41);
% 
% 
% for k= 1:params.num_cam
%   total_data_ct_single=[Train_mats.cam{k} Test_mats.cam{k}];
%   total_id_ct_multi=[Train_matm.camp{k} Test_matm.camp{k}];  
%   total_data_ct_multi=[Train_matm.cam{k} Test_matm.cam{k}];
%   
%   
%   [~,mod_multi_id_train]=find(total_id_ct_multi==train_rand_id');
%   [~,mod_multi_id_test]=find(total_id_ct_multi==test_rand_id');
% 
%   
%   
%   
%   Train_mats.camp{k}=train_rand_id;
%   Test_mats.camp{k}=test_rand_id;
%   
%   Train_mats.cam{k}=total_data_ct_single(:,random_index(1:25));
%   Test_mats.cam{k}=total_data_ct_single(:,random_index(26:41));
%   
%   Train_matm.camp{k}=total_id_ct_multi(mod_multi_id_train);
%   Test_matm.camp{k}=total_id_ct_multi(mod_multi_id_test);
%   
%   Train_matm.cam{k}=total_data_ct_multi(:,mod_multi_id_train);
%   Test_matm.cam{k}=total_data_ct_multi(:,mod_multi_id_test);
%   
%   
%   
%   
%   
%     
% end
% 
% end

%%




Metric_mat=[];

pairs=nchoosek(1:params.num_cam,2);


for i = 1:1 
params.ct_batch=i;
[Lmat_opt,stats]= opt_selpair(Train_mats,params);
stats.MA_mat=Lmat_opt;
if params.method ~=0
lb_gain_ratio = 100*lbg_ratio(Lmat_opt,params.master_budget);
else
lb_gain_ratio=0;
end
lb_pgain_ratio=stats.pos;



for j=1:size(pairs,1)
 full_ind_cam_a_train= Train_matm.camp{j,1};
 full_ind_cam_b_train= Train_matm.camp{j,2};
 single_ind_cam_a_train= Train_mats.camp{pairs(j,1)};
 single_ind_cam_b_train= Train_mats.camp{pairs(j,2)};
 
 full_ind_cam_a_test=Test_matm.camp{j,1};
 full_ind_cam_b_test=Test_matm.camp{j,2};
 Lmat_multi_ct_test=lmat_gen(full_ind_cam_a_test,full_ind_cam_b_test);

 augmented_ind_a_train=[];
 augmented_ind_b_train=[];
 Lmat_aug=[];
 

overlap_stat=overlap_stat_gen(stats.cam_ind);

[~,ct_ind_a_train]=find(full_ind_cam_a_train==(single_ind_cam_a_train(stats.indcam))');
[~,ct_ind_b_train]=find(full_ind_cam_b_train==(single_ind_cam_b_train(stats.indcam))');

augmented_ind_a_train=[augmented_ind_a_train ;   ct_ind_a_train];
augmented_ind_b_train=[augmented_ind_b_train  ;  ct_ind_b_train];


Lmat_multi_ct_train = replicate_lmat(Lmat_opt{j}, full_ind_cam_a_train(ct_ind_a_train)  ,full_ind_cam_b_train(ct_ind_b_train));
Lmat_aug=[Lmat_aug NaN(size(Lmat_aug,1),size(Lmat_multi_ct_train,2));NaN(size(Lmat_multi_ct_train,1),size(Lmat_aug,2)) Lmat_multi_ct_train ] ;




% temp_ux_a=[Train_matm.cam{pairs(j,1)}(:,ct_ind_a_train) Test_matm.cam{pairs(j,1)}];
% temp_ux_b=[Train_matm.cam{pairs(j,2)}(:,ct_ind_b_train) Test_matm.cam{pairs(j,2)}];
X_ct.train.data=[Train_matm.cam{j,1}(:,augmented_ind_a_train) Train_matm.cam{j,2}(:,augmented_ind_b_train)];
X_ct.train.Lmat=Lmat_aug;
X_ct.train.asmat=ones(size(X_ct.train.Lmat));
% X_ct.train.Lmat=lmat_gen(full_ind_cam_a_train,full_ind_cam_b_train);

X_ct.train.idxa=1:length(augmented_ind_a_train);
X_ct.train.idxb=length(X_ct.train.idxa)+1:size(X_ct.train.data,2);




X_ct.test.data=[ Test_matm.cam{j,1} Test_matm.cam{j,2}];
X_ct.test.Lmat=Lmat_multi_ct_test;
X_ct.test.asmat=ones(size(X_ct.test.Lmat));
X_ct.test.idxa=1:length(full_ind_cam_a_test);
X_ct.test.idxb=length(X_ct.test.idxa)+1:size(X_ct.test.data,2);
X_ct.test.inda=full_ind_cam_a_test;
X_ct.test.indb=full_ind_cam_b_test;









[~, Acc]= KISS_func(X_ct); 

ACC_mat{j,i}=Acc;
Req_vec(i)= 100*(stats.tot_req/stats.tot_var);

stats.overlap_stat=overlap_stat;
    
    
    
    
    
    
    
    
    
    
end

label_gain_vec(1,i)=lb_gain_ratio;
label_gain_vec(2,i)=lb_pgain_ratio;


end








end