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

clc; clear all; close all;
%DATA_OUT_DIR = fullfile('..','..','dataOut','cvpr','viper');
run('C:\Users\sroy\Google Drive\feat_\kissme\KISSME\toolbox\init.m');

%% Set up parameters

params.numCoeffs = 20000; %dimensionality reduction by PCA to 34 dimension
params.N = 20; %number of image pairs, 316 to train 316 to test
params.numFolds = 1; %number of random train/test splits
%params.saveDir = fullfile(DATA_OUT_DIR,'all');
params.pmetric = 0;

%% Load Features

%load(fullfile(DATA_OUT_DIR,'viper_features.mat'));
load mat_cam1set.mat
load mat_cam4set.mat

ux=[mat_cam1s(:,1:20) mat_cam4s(:,1:20)];
% [~,scores]=pca(ux','NumComponents',1000);
% ux=scores';
idxa=1:20;
idxb=21:40;

%% Cross-validate over a number of runs

% pair_metric_learn_algs = {...
%     LearnAlgoKISSME(params), ...
%     LearnAlgoMahal(), ...
%     LearnAlgoMLEuclidean(), ...
%     LearnAlgoITML(), ... 
%     LearnAlgoLDML(), ... 
%     LearnAlgoLMNN() ...  
%     };


pair_metric_learn_algs = {...
    LearnAlgoMLEuclidean(), ...
    };


[ ds ] = CrossValidateViper(struct(), pair_metric_learn_algs,ux(1:params.numCoeffs,:),idxa,idxb,params);

%% Plot Cumulative Matching Characteristic (CMC) Curves

names = fieldnames(ds);
for nameCounter=1:length(names)
   s = [ds.(names{nameCounter})];
   ms.(names{nameCounter}).cmc = cat(1,s.cmc)./(params.N/2);
   ms.(names{nameCounter}).roccolor = s(1).roccolor;
end

h = figure;
names = fieldnames(ms);
for nameCounter=1:length(names)
   hold on; plot(median(ms.(names{nameCounter}).cmc,1),'LineWidth',2, ...
       'Color',ms.(names{nameCounter}).roccolor);
end
  
title('Cumulative Matching Characteristic (CMC) Curves - VIPeR dataset');
box('on');
set(gca,'XTick',[0 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 150 200 250 300 350]);
ylabel('Matches');
xlabel('Rank');
ylim([0 1]);
hold off;
grid on;
legend(upper(names),'Location','SouthEast');

if isfield(params,'saveDir')
    exportAndCropFigure(h,'all_viper',params.saveDir);
    save(fullfile(params.saveDir,'all_data.mat'),'ds');
end
