clear all
close all
rand('state', 10);
addpath(genpath('UGM'));                                                            
addpath(genpath('libsvm-3.21/matlab'));
addpath(genpath('SFO'))

load cora.mat

NFeat = size(Feature,2)-2;  % No. of features
NData = size(Feature,1);    % No. of data points

% Divide the data in train-test
ind = crossvalind('kfold', Feature(:,end), 10);

trainset = Feature(ind~=1,:);
testset = Feature(ind==1,:);

[accuracy_vec, n_label_vec] = get_performance(trainset, testset, Link);
plot(cumsum(n_label_vec), accuracy_vec, '--o', 'LineWidth', 1);
xlabel('Number of manually labeled samples')
ylabel('Accuracy')