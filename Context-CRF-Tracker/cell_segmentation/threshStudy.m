% HIT F5 to run this file
%% Load the metric stored for analysis, this metric is used to find minimum
%  variance threshold
load threshStudyMetrics;
file_to_save_variables = './reuse_variables.mat';
[imgNum,threshNum,featNum] = size(threshMetrics);
%% set plotTrue = 0, get only values and no graph
% set plotTrue = 1, get values as well as graph
plotTrue=0;
xLab='Threshold';
yLab=['numCells';'avegArea';'varAreas'];
titleVec=['NumCells vs. Threshold';'AvegArea vs. Threshold';'Variance vs. Threshold'];
for l=1:imgNum,
    if(plotTrue==1)
        f(l)=figure;
        for k=1:3,
            subplot(1,featNum-1,k)
            plot(threshMetrics(l,1:threshNum,1), threshMetrics(l,1:threshNum,k+1));
            xlabel(xLab);
            ylabel(yLab(k,:));
            title(titleVec(k,:));
        end
    end
    
    b=threshMetrics(l,1:threshNum,featNum);
    minVar(l)=min(b);
    fVec=find(minVar(l)==b);
    if(length(fVec)>1)
        deltaThresh=mean(fVec)-1;  
        thresholdVec(l)=0.02+deltaThresh*0.005;
    else
        thresholdVec(l)=0.02+(fVec-1)*0.005;
    end
    l
end
%% Print the values
save(file_to_save_variables,'thresholdVec','-append');