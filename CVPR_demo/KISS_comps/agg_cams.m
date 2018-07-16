function [data ,id]= agg_cams(X)
data=[];
id=[];
for i=1:length(X.cam)
data=[data X.cam{i}];
id=[id X.camp{i}];
end

end