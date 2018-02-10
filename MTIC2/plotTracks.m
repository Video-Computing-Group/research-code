%% Code Author: Ahmed Tashrif Kamal - tashrifahmed@gmail.com
% http://www.ee.ucr.edu/~akamal
% no permission necessary for non-commercial use
% Date: 4/27/2013

%% 
function plotTracks(xa,mtic,icfgt,icfnn,jpdakcf,z,zId,zCount,Nc,Nt,t,T)
clf

hold on;
for i = 1:Nc
    if zCount(i)>0
        for n = 1:zCount(i)
            if zId{i}(n)>0 %true measurement
                h_obs = plot(z{i}(1,n),z{i}(2,n),'x','Color',[0,0,1],'MarkerSize',8);
            else %clutter
                h_clutter = plot(z{i}(1,n),z{i}(2,n),'.','Color',[0,0,0],'MarkerSize',8);
            end
        end
    end
end

i=1;

for j = 1:Nt
    plot(xa{j}(1,1:t),xa{j}(2,1:t),'g','LineWidth',2);
    plot(mtic(i,j).x(1,1:t),mtic(i,j).x(2,1:t),'r','LineWidth',1);
    %h_icfgt = plot(icfgt(i,j).x(1,1:t),icfgt(i,j).x(2,1:t),'m','LineWidth',1);
    %h_icfnn = plot(icfnn(i,j).x(1,1:t),icfnn(i,j).x(2,1:t),'b','LineWidth',1);
    %h_jpdakcf = plot(jpdakcf(i,j).x(1,1:t),jpdakcf(i,j).x(2,1:t),'k','LineWidth',1);
end

if t==T
    %false data for legend
    h_gt = plot([-100 -50],[-100 -50],'g','LineWidth',2);
    h_mtic = plot([-100 -50],[-100 -50],'r','LineWidth',1);
    %h_icfgt = plot([-100 -50],[-100 -50],'m','LineWidth',1);
    %h_icfnn = plot([-100 -50],[-100 -50],'b','LineWidth',1);    
    %h_jpdakcf = plot([-100 -50],[-100 -50],'k','LineWidth',1);
    h_obs = plot([-100 -50],[-100 -50],'x','Color',[0,0,1],'MarkerSize',8);
    h_clutter = plot([-100 -50],[-100 -50],'.','Color',[0,0,0],'MarkerSize',8);
    
    h_legend = legend([h_gt,h_mtic,h_obs,h_clutter],'GroundTruth','MTIC','Observation','Clutter');
    set(h_legend,'FontSize',12);
    
end




axis([0 500 0 500]);
axis square
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])

end