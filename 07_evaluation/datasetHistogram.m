function [] = datasetHistogram(data,binEdges,limitX,limitY)
    histogram(data,...
        'Normalization', 'probability','BinEdges',binEdges) %, 'NumBins', 10) %, 'FaceColor', colorCode)
    ylabel('Probability of occurrence')
    xlim(limitX)
    ylim(limitY)
    grid on

end