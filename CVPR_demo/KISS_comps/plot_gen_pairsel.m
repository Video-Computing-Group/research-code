function [] =plot_gen_pairsel(struct,num_pairs)
close all
methods=struct.methods;
method_marker=struct.marker;
for i=num_pairs: num_pairs
    figure;
    for j= 1: length(methods)-1
    if(length(struct.acc{1,j})<1)
        continue;
    end
    plot((struct.req{1,j})',100*(struct.acc{1,j}(:,i)),method_marker{j},'LineWidth',2);

    hold on
        
        
        
    end
    plot([ 0 100],100*[struct.full_set.acc(i) struct.full_set.acc(i)],method_marker{j+1},'LineWidth',2);
    legend(methods{1},methods{2},methods{3},methods{4})
    axis([0 100 0 100])
    xlabel('Percentage of Manual Labeling')
    ylabel('Accuracy')
    hold off 

    
    
end



end