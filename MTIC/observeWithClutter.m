%% Code Author: Ahmed Tashrif Kamal - tashrifahmed@gmail.com
% http://www.ee.ucr.edu/~akamal
% no permission necessary for non-commercial use
% Date: 4/27/2013

%% 
function [zt,zIdt,zCountt]=observeWithClutter(lamda,Nc,Nt,xa,H,R,FOV,T)
zt = cell(T,1);
zIdt = cell(T,1);
zCountt = cell(T,1);

meanc = [0;0];

for t=1:T
    
    z = cell(Nc,1);
    zId = cell(Nc,1);
    zCount = zeros(Nc,1);
    
    for i = 1:Nc
        z{i}=[];
        zId{i}=[];
        for j=1:Nt
            
            % generate observations
            ztemp = H*xa{j}(:,t);
            if pointInRect(ztemp,FOV(:,i))
                z{i} = [z{i}, ztemp+mvnrnd(meanc,R)'];
                zId{i} = [zId{i},j];
                zCount(i) = zCount(i) + 1;
            end
            
            %add clutter
            sensorPos = FOV(1:2,i);
            sensorSize = [FOV(3,i)-FOV(1,i);FOV(4,i)-FOV(2,i)];
            numClutters = poissrnd(lamda);
            for c = 1:numClutters
                z{i} = [z{i}, sensorSize.*rand([2,1])+sensorPos     ];
                zId{i} = [zId{i},0];
                zCount(i) = zCount(i) + 1;
            end
            
        end
    end
    zt{t} = z;
    zIdt{t} = zId;
    zCountt{t} = zCount;
end

end