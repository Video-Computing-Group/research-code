clear 
clc
params.num_cam=3;
params.num_person=35;
params.non_overlap=1;
params.non_overlap_lim=35;
params.method=0;
params.budget=10;
params.randomize_split=1;
params.reg_spec=8;
num_fold_val=1;
master_budget_t=15; %%% Change this to change labeling budget

rng('default')
params.itcl=450;
ACC_mat1_store=cell(nchoosek(params.num_cam,2),1);
Req_vec1_store=0;
gain_vec1_store=0;
for k=1:1
params.num_folds=1;
[ACC_mat1,Req_vec1,gain_vec1,stats] = Ward_sim(params);
%ACC_mat1_store=ACC_mat1_store+ACC_mat1';
ACC_mat1_store=cell_add_1d(ACC_mat1_store,ACC_mat1);
Req_vec1_store=Req_vec1+Req_vec1_store;
gain_vec1_store=gain_vec1+gain_vec1_store;


end
%ACC_mat1_avg=ACC_mat1_store/k;
%num_folds=5;
ACC_mat1_avg=cell_normalizer_1d(ACC_mat1_store,k);
Req_vec1_avg=Req_vec1_store/k;


 budget_mat=0.7*master_budget_t;

% budget_mat=[18];
master_budget=(stats.tot_var/100)*[master_budget_t];


%%% exact solution
for method_iter=1:4
params.method=method_iter;

Method_acc_mat=[];
Method_req_mat=[];
Method_gain_mat=[];
params.budget=budget_mat;
params.master_budget=floor(master_budget);
% if method_iter==3
%     params.full_budget=round(stats.tot_var*(Method.req{1,1}(j))/100);
%     
% end
params.full_budget=master_budget;
if method_iter==4
    params.full_budget=(stats.tot_var/100)*master_budget_t;
end

ACC_mat_store=cell( nchoosek(params.num_cam,2),1);
Req_vec_store=0;
gain_vec_store=0;

rng('default')
if method_iter==5
    num_folds=num_fold_val;
else
   num_folds=1;     
end
for i=1:num_folds
params.num_folds=1;

[ACC_mat2,Req_vec2,gain_vec2,stats] = Ward_sim(params);
ACC_mat_store=cell_add_1d(ACC_mat_store,ACC_mat2');
Req_vec_store=Req_vec2+Req_vec_store;
gain_vec_store=gain_vec2+gain_vec_store;

Method.FL{method_iter,i}=stats.MA_mat;
Method.LC{method_iter,i}=stats.LC;
end
% ACC_mat_Avg=ACC_mat_store/i;
% Req_vec_Avg=Req_vec_store/i;
ACC_mat_Avg=cell_normalizer_1d(ACC_mat_store,i);
Req_vec_Avg=Req_vec_store/i;
gain_vec_Avg=gain_vec_store/i;

Method_acc_mat=[Method_acc_mat; ACC_mat_Avg];
Method_req_mat=[Method_req_mat; Req_vec_Avg];
Method_gain_mat=[Method_gain_mat; gain_vec_Avg];




Method.acc{method_iter}=Method_acc_mat;
Method.req{method_iter}=Method_req_mat;
Method.gain{method_iter}=Method_gain_mat;





end

Method.full_set.acc=ACC_mat1_avg;
Method.full_set.req=Req_vec1_avg;
Method.marker{1}='-rs';
Method.marker{2}='--cs';
Method.marker{3}='--b*';
Method.marker{4}='-m';



%Method.acc{1}=ACC_mat1_avg;

Method.methods{1}='Exact';
Method.methods{2}='Greedy';
Method.methods{3}='Half_apx';
Method.methods{4}='Baseline';


Method.parameters.mb=master_budget_t;
Method.parameters.sb=budget_mat;
Method.parameters.reg=params.reg_spec;
Method.parameters.folds=num_folds;
Method.parameters.sp=params.non_overlap_lim;

Gain_res=[Method.gain{1} Method.gain{2} Method.gain{3} Method.gain{4}];
rmtrc{1}='Total Labels(%)';
rmtrc{2}='Positive Labels (%)';
clc
Label_Recovered = array2table(Gain_res,'variableNames',Method.methods,'RowNames',rmtrc)
ot=method_acc_avg_cam(Method);

Method.methods{5}='Fullset';

accl2{1}='Rank-1 accuray(%)';
Rank_1_accuracy = array2table(ot,'variableNames',Method.methods,'RowNames',accl2)
FL_store=Method.FL;
MA_store=Method.LC;
load Ward_dataset_folds3.mat
Data_set= Ward_dataset_folds2;
clearvars -except MA_store FL_store Data_set
% 
% strsv=strcat('Rn_',num2str(params.non_overlap_lim),'_',num2str(master_budget_t),'_',num2str(round(budget_mat)),'.mat');
% save(strsv,'Method')
