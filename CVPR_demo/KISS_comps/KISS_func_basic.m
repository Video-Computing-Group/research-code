%% Person Re-Identification on the VIPeR dataset [6]
%
% [6] D. Gray, S. Brennan, and H. Tao. Evaluating appearance
% models for recongnition, reacquisition and tracking. In Proc.
% IEEE Intern.Workshop on Performance Evaluation of Tracking
% and Surveillance, 2007.
%
% See also http://vision.soe.ucsc.edu/?q=node/178
%
% Features:
%
% HSV, Lab histograms and LBPs [16] to describe color and texture of the
% overlapping blocks (size 8x16 and stride of 8x8 pixels). The image 
% descriptor is a concatenation of the local ones. Using PCA the descriptor
% is projected onto a 34 dimensional subspace.
% 

function [Metric_mat, Acc_vec]= KISS_func_basic(X)

% clc; clear ; close all;
% load mat_cam1set.mat
% load mat_cam2set.mat
% load mat_cam3set.mat
% load mat_cam4set.mat
% load sel_return
% 
% 
% 
% 
% feat{1}=mat_cam1s;
% feat{2}=mat_cam2s;
% feat{3}=mat_cam3s;
% feat{4}=mat_cam4s;
% 
% 
% pairs=nchoosek([1 2 3 4],2);
%%%

params.numCoeffs =100; %dimensionality reduction by PCA to 34 dimension
%params.N = 41; %number of image pairs, 316 to train 316 to test
params.numFolds = 1; %number of random train/test splits
params.pmetric = 0;
params.pca=0;
params.less_neg=1;
%%%

for i=3:3
% ct_cama=feat{pairs(i,1)};  
% ct_camb=feat{pairs(i,2)};  

% asmat_sel=lmat_structc2{i};
% asmat=ones(size(lmat_structc2{i}));

% params.idx{1}=cam_ind{pairs(i,1)};  
% params.idx{2}=cam_ind{pairs(i,2)};  
% params.asmat=asmat;

% params.lmat=(tlmatstruct{i}).*asmat_sel;
% idxa=1:41;
% idxb=42:82;
% idxa=params.idxa;
% idxb=params.idxb;

%%%%
pair_metric_learn_algs = {...
     LearnAlgoKISSME(), ...
     };
%  pair_metric_learn_algs = {...
%       LearnAlgoMahal(), ...
%      };
%%%%
% ux=[ct_cama ct_camb];
% [ ds ] = CrossValidate_multi(struct(), pair_metric_learn_algs,ux(1:end,:),idxa,idxb,params);
 [ ds ] = CrossValidate_multi_basic(struct(), pair_metric_learn_algs,X,params);

    
     
    
end

 names = fieldnames(ds);
% for nameCounter=1:length(names)
%    s = [ds.(names{nameCounter})];
%    ms.(names{nameCounter}).cmc = cat(1,s.cmc)./16;
%    ms.(names{nameCounter}).roccolor = s(1).roccolor;
% end
% acc=mean(ms.kissme.cmc(:,1))
% 
% 
% 
% 
% 
Metric_mat=[];
Acc_vec=ds.accuracy;
 end
