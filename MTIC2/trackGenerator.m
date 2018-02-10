%% Code Author: Ahmed Tashrif Kamal - tashrifahmed@gmail.com
% http://www.ee.ucr.edu/~akamal
% no permission necessary for non-commercial use
% Date: 4/27/2013

%% Initialization
disp('Generating Ground Truth Tracks...');
xa = cell(Nt,1); %groundtruth state
j = 1;

while(j<=Nt)
    flag_restart = 0;
    
    %% initial ground truth position
    speed = mag_min+(mag_max-mag_min)*rand;
    direction = angle_min+(angle_max-angle_min)*rand;
    x = gridxMin + (gridxMax-gridxMin) * rand;
    y = gridyMin + (gridyMax-gridyMin) * rand;
    xa{j} = [x ; y ; speed*cos(direction) ; speed*sin(direction)  ];
    
    for t = 1:T
        %% generate ground truth tracks
        if t>1
            xa{j}(:,t) = Phi*xa{j}(:,t-1)+mvnrnd(zeros(p,1),Q)';
            
            % keep iterating only if within the bounding box
            if xa{j}(1,t)> 500 || xa{j}(1,t) < 0 || xa{j}(2,t) > 500 || xa{j}(2,t) < 0
                flag_restart = 1;
                break;
            end
        end
    end %timestep t
    
    if flag_restart
        continue;
    elseif norm( xa{j}(1:2,1) - xa{j}(1:2,T) ) < 200 %discard small tracks
        continue;
    else
        j = j + 1;
    end
    
    
    
end %end j

save ('groundTruthTracks','xa');
