function [ result ] = calcMCMC_multi( M, data, idxa, idxb, idxtest,tlabel )
%%% multi-shot CMC curve generator 


% dist = sqdist(data(:,idxa(idxtest)), data(:,idxb(idxtest)),M);
% 
% result = zeros(1,size(dist,2));
% for pairCounter=1:size(dist,2)
%     distPair = dist(pairCounter,:);  
%     [tmp,idx] = sort(distPair,'ascend');
%     result(idx==pairCounter) = result(idx==pairCounter) + 1;
% end
% 
% tmp = 0;
% for counter=1:length(result)
%     result(counter) = result(counter) + tmp;
%     tmp = result(counter);
% end

%%

dist=Distmat((data(:,idxa(idxtest)))', (data(:,idxb(idxtest)))',M);
[~,inds]=sort(dist,2,'ascend');
result=zeros(1,length(idxtest));
for i=1:size(dist,1)
predl=tlabel(inds(i,:));
indct=find(predl==tlabel(i),1);
result(1,indct:end)=result(1,indct:end)+1;
    
    
end



end

