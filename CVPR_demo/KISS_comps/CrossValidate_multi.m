function [ ds ] = CrossValidate_multi(ds, learn_algs, X, params)
%function [ds,runs]=CrossValidateViper(ds,learn_algs,X,idxa,idxb,params)
%   
% Input:
%
%  ds - data struct that stores the result
%  learn_algs - algorithms that are used for cross validation
%  X - input matrix, each column is an input vector [DxN*2]. N is the
%  number of pairs.
%  idxa - index of image A in X [1xN]
%  idxb - index of image B in X [1xN]
%  params.N - number of pairs
%  params.numFolds - number of runs
%
% Output:
%
%  ds - struct [1xnumFolds] that contains the result
%  runs - struct [1xnumFolds] that contains the train test split
%  runs(c).perm - random permutation of run c
%  runs(c).idxtrain - train index
%  runs(c).idxtest - test index
%
% See also CrossValidatePairs
%
% copyright by Martin Koestinger (2011)
% Graz University of Technology
% contact koestinger@icg.tugraz.at
%
% For more information, see <a href="matlab: 
% web('http://lrs.icg.tugraz.at/members/koestinger')">the ICG Web site</a>.

for c=1:1
    %clc; fprintf('VIPeR run %d of %d\n',c,params.numFolds);
idxa=X.train.idxa;
idxb= X.train.idxb;


    % split in equal-sized train and test sets
%     idxtrain1 = params.idx{1};
%     idxtrain2 = params.idx{2};
%     asmat=params.asmat;
%     lmat=params.lmat;
    pca_flag=params.pca;
  %%  tlabel=params.tlabel;
    
    %%
%     idxtest  = 26:41;
%     idxtrain=1:25;

    if pca_flag==1
    X_tr=X.train.data;%X(:,[idxa(idxtrain1) idxb(idxtrain2)]);
    X_te=X.test.data;%X(:,[idxa(idxtest) idxb(idxtest)]);
    
    XX_tr=zeros(params.numCoeffs,size(X_tr,2));%zeros(params.numCoeffs,size(X,2));
    XX_te=zeros(params.numCoeffs,size(X_te,2));
    
    
    
    [COEFF,SCORE]=pca(X_tr','NumComponents',params.numCoeffs);
    SCORE2=bsxfun(@minus,X_te',mean(X_tr'))*COEFF;
    %XX1(:,[idxa(idxtrain1) idxb(idxtrain2)])=SCORE';
    %XX1(:,[idxa(idxtest) idxb(idxtest)])=SCORE2';
    %X=XX1;
    %clear XX1
    XX_tr=SCORE';
    XX_te=SCORE2';
    
    else
    XX_tr=X.train.data;
    XX_te=X.test.data;
    end
   %%

    [Idx_A ,Idx_B,La1]= completeind2ff(X.train.asmat,X.train.Lmat,1:length(idxa),1:length(idxb),idxa,idxb);
%     if params.less_neg==1
%         neg_lim=600;
%         Idx_A=Idx_A(1:neg_lim);
%         Idx_B=Idx_B(1:neg_lim);
%         La1=La1(1:neg_lim);
% 
%         
%     end
% lim_pairs_kiss=200000;
% if length(La1)>lim_pairs_kiss
%   Idx_A=Idx_A(1:lim_pairs_kiss);
%   Idx_B=Idx_B(1:lim_pairs_kiss);
%   La1=La1(1:lim_pairs_kiss);
% 
% end
   if(length(find(La1==1))>1) 
    % train on first half
    for aC=1:length(learn_algs)
        cHandle = learn_algs{aC};
        fprintf('    training %s ',upper(cHandle.type));
         s = learnPairwise(cHandle,XX_tr,Idx_A,Idx_B,(La1));
        if ~isempty(fieldnames(s))
            fprintf('... done in %.4fs\n',s.t);
            ds(c).(cHandle.type) = s;
        else
            fprintf('... not available');
        end
    end
       
    % test on second half
    names = fieldnames(ds(c));
    for nameCounter=1:length(names)       
        fprintf('    evaluating %s ',upper(names{nameCounter}));
        %ds(c).(names{nameCounter}).cmc = calcMCMC(ds(c).(names{nameCounter}).M, XX_te,X.test.idxa,X.test.idxb,1:16);
        [~,acc]= nn1clf_k((XX_te(:,X.test.idxa)),(XX_te(:,X.test.idxb)),X.test.inda,X.test.indb,ds(c).(names{nameCounter}).M);
       %acc(1)
        fprintf('... done \n');
        ds.accuracy=acc;
    end
   else
        ds.accuracy=0
   
 end

end

end