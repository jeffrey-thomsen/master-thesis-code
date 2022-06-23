% plot histograms and heat maps for showcasing target and interferer
% precision and recall rates
function [] = plotDOAHistHeat(precisionTarget, precisionInterf, recallTarget, recallInterf)
    
%     figure('Position',[50 50 2.65*576 2.3*432]);
    figure('Position',[50 50 1.325*576 1.15*432]);
    subplot = @(m,n,p) subtightplot(m,n,p,[0.1 0.1], [0.15 0.03], [0.1 0.12]);

    subplot(1,2,1)
    datasetHistogram(precisionTarget,0:0.025:1,[0 1],[0 0.9])
    xlabel('Rate (1)')
    grid on
    set(gca,'FontSize',14);
%     title('Target precision');
    
    subplot(1,2,2)
    datasetHeatmap(precisionTarget,'',[]);
    %     title('Target precision');
    set(gca,'FontSize',14);
    
    figure('Position',[50 50 1.325*576 1.15*432]);
    subplot(1,2,1)
    datasetHistogram(precisionInterf,0:0.025:1,[0 1],[0 0.9])
    xlabel('Rate (1)')
    grid on
%     title('Interfer precision');
    set(gca,'FontSize',14);
    
    subplot(1,2,2)
    datasetHeatmap(precisionInterf,'',[]);
    %     title('Interfer precision');
    set(gca,'FontSize',14);
    
    figure('Position',[50 50 1.325*576 1.15*432]);
    subplot(1,2,1)
    datasetHistogram(recallTarget,0:0.025:1,[0 1],[0 0.9])
    xlabel('Rate (1)')
    grid on
%     title('Target recall');
    set(gca,'FontSize',14);
    
    subplot(1,2,2)
    datasetHeatmap(recallTarget,'',[]);
    %     title('Target recall');
    set(gca,'FontSize',14);
    
    figure('Position',[50 50 1.325*576 1.15*432]);
    subplot(1,2,1)
    datasetHistogram(recallInterf,0:0.025:1,[0 1],[0 0.9])
    xlabel('Rate (1)')
    grid on
%     title('Interfer recall');
    
    set(gca,'FontSize',16);
    subplot(1,2,2)
    datasetHeatmap(recallInterf,'',[]);
    %     title('Interfer recall');
    set(gca,'FontSize',16);

end