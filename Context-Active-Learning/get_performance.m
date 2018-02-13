function [accuracy_vec, n_label_vec] = get_performance(trainset, testset, link)

% The training set is divided into batches of data and the first batch of
% data is considered to be fully manually labeled. The batches thereafter
% are considered to be unlabeled and a subset of each of those batches are
% selected for manual labeling and added to the labeled set and then the
% current model is updated

n_batch = 10;
accuracy_vec = zeros(1,n_batch);
n_label_vec = zeros(1,n_batch);
idx = crossvalind('kfold', trainset(:,end), n_batch);
categories = unique([trainset(:,end);testset(:,end)]);
edgepot = zeros(length(categories));

batchdata = trainset(idx==1,:);
labeleddata = batchdata;
n_label_vec(1) = size(batchdata,1);

model = svmtrain(batchdata(:,end), batchdata(:,2:end-1), '-t 0 -c 1 -b 1 -q');
disp('Successfully created initial model !!')
[~, accfirst, ~] = svmpredict(testset(:,end), testset(:,2:end-1), model, '-b 1 -q');
accuracy_vec(1) = accfirst(1);

edgepot = get_edgepot(edgepot, link, batchdata(:,1), batchdata(:,end));

% Iterate over number of batches of incoming data
for i = 2:n_batch
    disp(['Processing batch ' num2str(i) ' of ' num2str(n_batch)])
    batchdata = trainset(idx==i,:);
    [labeleddata, acc, n_label, newedgepot, newmodel] = get_batchperformance(labeleddata, batchdata, testset, link, edgepot, model);
    accuracy_vec(i) = acc(1);
    n_label_vec(i) = n_label;
    edgepot = newedgepot;
    model = newmodel;
end

