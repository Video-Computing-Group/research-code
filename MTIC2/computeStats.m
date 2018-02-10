%% Code Author: Ahmed Tashrif Kamal - tashrifahmed@gmail.com
% http://www.ee.ucr.edu/~akamal
% no permission necessary for non-commercial use
% Date: 4/27/2013

%% 

function [ME_mtic,ME_icfgt,ME_icfnn,ME_jpdakcf,SDE_mtic,SDE_icfgt,SDE_icfnn,SDE_jpdakcf,e_mtic,e_icfgt,e_icfnn,e_jpdakcf]...
    = computeStats(xa,mtic,icfgt,icfnn,jpdakcf)
T = size(mtic(1,1).x,2);
[Nc,Nt] = size(mtic);
N = Nc*Nt*T;

%% individual errors
e_mtic = zeros(Nc,Nt,T);
e_icfgt = zeros(Nc,Nt,T);
e_icfnn = zeros(Nc,Nt,T);
e_jpdakcf = zeros(Nc,Nt,T);

for i=1:Nc
    for j=1:Nt
        for t=1:T
            e_mtic(i,j,t) = norm(xa{j}(1:2,t)-mtic(i,j).x(1:2,t));
            e_icfgt(i,j,t) = norm(xa{j}(1:2,t)-icfgt(i,j).x(1:2,t));
            e_icfnn(i,j,t) = norm(xa{j}(1:2,t)-icfnn(i,j).x(1:2,t));
            e_jpdakcf(i,j,t) = norm(xa{j}(1:2,t)-jpdakcf(i,j).x(1:2,t));
        end
    end
end


%% mean error
ME_mtic = sum(e_mtic(:))/N;
ME_icfgt = sum(e_icfgt(:))/N;
ME_icfnn = sum(e_icfnn(:))/N;
ME_jpdakcf = sum(e_jpdakcf(:))/N;
 

%% sqrerror
esqr_mtic = e_mtic.*e_mtic;
esqr_icfgt = e_icfgt.*e_icfgt;
esqr_icfnn = e_icfnn.*e_icfnn;
esqr_jpdakcf = e_jpdakcf.*e_jpdakcf;


%% standard deviation of error
SDE_mtic = sqrt( sum(esqr_mtic(:)) /N - ME_mtic*ME_mtic);
SDE_icfgt = sqrt( sum(esqr_icfgt(:)) /N - ME_icfgt*ME_icfgt);
SDE_icfnn = sqrt( sum(esqr_icfnn(:)) /N - ME_icfnn*ME_icfnn);
SDE_jpdakcf = sqrt( sum(esqr_jpdakcf(:)) /N - ME_jpdakcf*ME_jpdakcf);
 
disp(['MTIC: MeanError: ' num2str(ME_mtic) ', Standard Deviation: ' num2str(SDE_mtic) ]);
disp(['ICF-GT: MeanError: ' num2str(ME_icfgt) ', Standard Deviation: ' num2str(SDE_icfgt) ]);
disp(['ICG-NN: MeanError: ' num2str(ME_icfnn) ', Standard Deviation: ' num2str(SDE_icfnn) ]);
disp(['JPDA-KCF: MeanError: ' num2str(ME_jpdakcf) ', Standard Deviation: ' num2str(SDE_jpdakcf) ]);

save('MTIC_Results','ME_mtic','ME_icfgt','ME_icfnn','ME_jpdakcf','SDE_mtic','SDE_icfgt','SDE_icfnn',...
    'SDE_jpdakcf','e_mtic','e_icfgt','e_icfnn','e_jpdakcf');

end