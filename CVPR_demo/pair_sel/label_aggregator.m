function [out1, out2]=label_aggregator(tempx1,tempx2,pair_1label,pair_2label)
% pair_2labelt(pair_2labelt==0)=-1;
% A1=pair_2label.*pair_2labelt;
out1=tempx1*pair_2label;

% pair_1labelt(pair_1labelt==0)=-1;
% A3=pair_1label.*pair_1labelt;
out2=pair_1label*tempx2;









end