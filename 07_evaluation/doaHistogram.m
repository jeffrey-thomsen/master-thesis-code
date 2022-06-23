function [] = doaHistogram(varargin)

    azDeg = varargin{1};
    if nargin == 2
        colorCode = varargin{2};
    elseif nargin == 1
        colorCode = '#D95319';
    else
        error('Must specify 1 or 2 input arguments.')
    end

    histogram(nonzeros(azDeg(abs(azDeg)<90)),...
        'Normalization', 'probability', 'NumBins', 90, 'FaceColor', colorCode)
%     title(['DOA estimation histogram'])%,['total no. of subband samples: ',num2str(length(testSignal)*nBands)])
%     ylabel({'No. of occurences';'in subband samples'})
%     ylabel('Probability of occurrence')
%     xlabel('$\hat{\varphi} (^{\circ})$','Interpreter','latex')
    xticks([-90 -60 -30 0 30 60 90])
%     ylim([0 0.25])
    grid on
    set(gca,'FontSize',12);

end