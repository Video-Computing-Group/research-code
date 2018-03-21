%%-------------------------------------------------------------------------
% Main program to find the summaries (key shots) in Collaborative Video Summarization 
% Input: multiple video feature matrices (shot level C3D features) - one is
% for the target video and other is the side information that will be used
% to assist summarizing the target video
% Output: representative list of shots
%--------------------------------------------------------------------------
%%-------------------------------------------------------------------------
clc; clear all; close all
%% Load the feature matrices
tic;
fprintf('Loading the video feature matrices:\n');

% To make it simple, lets assume x and y are two feature matrices 
% x - target video to be summarized
% y - other topic related videos - combined to form a single feature matrix

load('Y:\Rameswar\ECCV_2016\Features\CoSum\C3DFeatures\BJ\C3D_F1.mat'); % Target Video
X = feature_mat;
load('Y:\Rameswar\ECCV_2016\Features\CoSum\C3DFeatures\BJ\C3D_F2.mat'); % Auxilliary Videos - concatenate all to a single video
Y = feature_mat;

nvideo = 2;
cellarray = cell(1,nvideo);
cellarray{1} = X;
cellarray{2} = Y;

% Feature Zero mean
for i=1:nvideo
    cellarray{i} = cellarray{i} - repmat(mean(cellarray{i},2),1,size(cellarray{i},2));
end

fprintf('Successfully Completed !!!\n');
fprintf('-----------------------------------\n');

%% Computing the similarities for representation diversity regularizations
fprintf('Computing the pairwise video similarities:\n');
simcellarray = cell(nvideo,nvideo);

% Use cosine similarity, others can also be used
simtype = 'Please enter the type of similarity to be used: 1-cosine, 2-pearson cor coeff, 3-gaussian+SLH?';
sim = input(simtype);

if sim == 1
% Computing Cosine similarity

for k1=1:nvideo
    for k2=1:nvideo
        D = zeros(size(cellarray{k1},2),size(cellarray{k2},2));
        for k3=1:size(cellarray{k1},2)
            for k4=1:size(cellarray{k2},2)
               D(k3,k4) = dot(cellarray{k1}(:,k3),cellarray{k2}(:,k4))/(norm(cellarray{k1}(:,k3))*norm(cellarray{k2}(:,k4)));
               simcellarray{k1,k2}= max(D,0); % Restricting to nonnegative
               %simcellarray{k1,k2}= D;
            end
        end
        clear D;
    end
end

% Computing the Pearson's correlation coefficient

elseif sim == 2
for k1=1:nvideo
    for k2=1:nvideo
        D = zeros(size(cellarray{k1},2),size(cellarray{k2},2));
        for k3=1:size(cellarray{k1},2)
            for k4=1:size(cellarray{k2},2)
               co = cov(cellarray{k1}(:,k3),cellarray{k2}(:,k4));
               D(k3,k4) = co(2)/(std(cellarray{k1}(:,k3))*std(cellarray{k2}(:,k4)));
               simcellarray{k1,k2}= max(D,0); % Restricting to nonnegative
               %simcellarray{k1,k2}= D;
            end
        end
        clear D;
    end
end

% Computing similarity via Gaussian similarities and SLH algorithm

else
for k1=1:nvideo
    for k2=1:nvideo
        Dis = pdist2(transpose(cellarray{k1}),transpose(cellarray{k2}));
        sigma = 0.2*max(max(Dis));
        Dg = exp(-Dis.^2 ./ (2*sigma^2));
        %Applying the SLH
        [U,S,V] = svd(Dg);
        S(logical(eye(size(S)))) = 1;
        D = U*S*V';
        simcellarray{k1,k2}= max(D,0); % Restricting to nonnegative
        %simcellarray{k1,k2}= max(Dg,0);
        %simcellarray{k1,k2}= D;
        clear D;
    end
end
end

% To make it simple, lets assume D1 and D2 are the two similarity matrices
% D1 - weight matrix measuring the pair-wise similarity of shots in X
% D2 - weight matrix measuring the pair-wise similarity of shots in X and Y

D1 = simcellarray{1,1};
D2 = simcellarray{1,2};

fprintf('Successfully Completed !!!\n');
fprintf('-----------------------------------\n');

%% Half-quadratic optimization Parameters
fprintf('Intializing Optimization Parameters: \n');

q = 2;
%alpha = 0.5; % typically alpha in [0,1]
alphainput = 'Please enter the value of regularization parameter - alpha?';
alpha = input(alphainput);

%lambda_d = 0.01; % typically lambda_d in [0.1, 0.01, 0.001, 0.0001]
lambdadinput = 'Please enter the value of regularization parameter - lambda_d?';
lambda_d = input(lambdadinput);

%gamma = 50; % typically gamma in [10,100]
gammainput = 'Please enter the value of regularization parameter - gamma?';
gamma = input(gammainput);

lambda_s = GetLambdaData(X)*1/gamma;
beta = lambda_s;

maxIter = 500;

%% Half-quadratic optimization
% Ref: Efficient and Robust Feature Selection via Joint L2,1-Norms
% Minimization, NIPS'10

fprintf('Start of Optimization: \n');

Z_c = HQCvs(cellarray{1},cellarray{2},D1,D2,alpha,lambda_s,lambda_d,beta,maxIter);


%% Finding the representative shots from the sparse coefficient matrix

fprintf('Finding the representative shots:\n');

t1 = 0.90; 
rank_list = findrank(Z_c,t1);

t2 = 0.95;
replist = rmRep(rank_list,cellarray{1},t2); % Ref: See All by Looking at A Few, CVPR'12.

%% 

fprintf('Select Representative Shots from replist: \n');

fprintf('Successfully Completed !!!\n');
fprintf('-----------------------------------\n');


