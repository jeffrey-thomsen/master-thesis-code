function [] = plotDOAHistHeat(precisionTarget, precisionInterf, recallTarget, recallInterf)
    
    figure('Position',[10 10 2.65*576 2.3*432]);
    
    subplot(2,4,1)
    datasetHistogram(precisionTarget,0:0.025:1,[0 1],[0 0.9])
    xlabel('Rate (1)')
    grid on
    set(gca,'FontSize',11);
    title('Target precision');
    
    subplot(2,4,2)
    datasetHeatmap(precisionTarget,'',[]);
    %     title('Target precision');
    set(gca,'FontSize',11);
    
    subplot(2,4,3)
    datasetHistogram(precisionInterf,0:0.025:1,[0 1],[0 0.9])
    xlabel('Rate (1)')
    grid on
    title('Interfer precision');
    set(gca,'FontSize',11);
    
    subplot(2,4,4)
    datasetHeatmap(precisionInterf,'',[]);
    %     title('Interfer precision');
    set(gca,'FontSize',11);
    
    subplot(2,4,5)
    datasetHistogram(recallTarget,0:0.025:1,[0 1],[0 0.9])
    xlabel('Rate (1)')
    grid on
    title('Target recall');
    set(gca,'FontSize',11);
    
    subplot(2,4,6)
    datasetHeatmap(recallTarget,'',[]);
    %     title('Target recall');
    set(gca,'FontSize',11);
    
    subplot(2,4,7)
    datasetHistogram(recallInterf,0:0.025:1,[0 1],[0 0.9])
    xlabel('Rate (1)')
    grid on
    title('Interfer recall');
    
    set(gca,'FontSize',11);
    subplot(2,4,8)
    datasetHeatmap(recallInterf,'',[]);
    %     title('Interfer recall');
    set(gca,'FontSize',11);

end