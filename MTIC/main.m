%% Code Author: Ahmed Tashrif Kamal - tashrifahmed@gmail.com
% http://www.ee.ucr.edu/~akamal
% no permission necessary for non-commercial use
% Date: 4/27/2013

%% 
clc;
clear all;
%clear classes;

disp('start');
generateFreshData = false;

if generateFreshData
    run LoadParameters
    run trackGenerator
    run sensorGenerator
    run PriorErrorGenerator
    run observationGenerator
else
    load savedParameters.mat
    load groundTruthTracks.mat
    load groundTruthTracks.mat
    load FOV.mat
    load priorError.mat
    load savedObservations.mat
end

%initialize plot
drawSensorNetwork(FOV,E);
h_fig_track=figure(2);
clf

%initialize filters
mtic = MTIC(Nc,Nt,p,T,eta,xa,Pinv);
icfgt = ICF_GT(Nc,Nt,p,T,eta,xa,Pinv);
icfnn = ICF_NN(Nc,Nt,p,T,eta,xa,Pinv);
jpdakcf = JPDA_KCF(Nc,Nt,p,m,T,eta,xa,P);


for t=1:T %timestep
    
    %% get measurements
    z = zt{t};
    zId = zIdt{t};
    zCount = zCountt{t};
    
    %% MTIC
    mtic.dataAssociation(Nc,Nt,p,m,z,zCount,R,H,FOV,Pg,lamdaf);
    mtic.prepData(Nc,Nt);
    mtic.consensus(Nc,Nt,p,K,eps,E);
    mtic.estimate(Nc,Nt,p,t,Pinv,Phi,Q);
    
    %% ICF-GT
    icfgt.dataAssociation(Nc,Nt,p,z,zId,zCount,Rinv,H);
    icfgt.prepData(Nc,Nt);
    icfgt.consensus(Nc,Nt,p,K,eps,E);
    icfgt.estimate(Nc,Nt,p,t,Phi,Q);
    
    %% ICF-NN
    icfnn.dataAssociation(Nc,Nt,p,z,zCount,Rinv,H,FOV);
    icfnn.prepData(Nc,Nt);
    icfnn.consensus(Nc,Nt,p,K,eps,E);
    icfnn.estimate(Nc,Nt,p,t,Phi,Q);
    
    %% JPDA-KCF
    jpdakcf.dataAssociation(Nc,Nt,p,m,z,zCount,R,H,FOV,Pg,lamdaf);
    jpdakcf.commOnce(Nc,Nt,p,E);
    jpdakcf.estimate(Nc,Nt,p,t,eps,E);
    jpdakcf.consensus(Nc,Nt,p,K,E,t);
    jpdakcf.predict(Nc,Nt,t,P,Phi,Q);
    
    %% Plot Tracking
    figure(h_fig_track);
    plotTracks(xa,mtic,icfgt,icfnn,jpdakcf,z,zId,zCount,Nc,Nt,t,T);
    pause(.001);
end %t

[ME_mtic,ME_icfgt,ME_icfnn,ME_jpdakcf,SDE_mtic,SDE_icfgt,SDE_icfnn,SDE_jpdakcf,e_mtic,e_icfgt,e_icfnn,e_jpdakcf]...
    = computeStats(xa,mtic,icfgt,icfnn,jpdakcf);

% save plot
box on;
saveas(h_fig_track,'MTIC_Results.pdf');

