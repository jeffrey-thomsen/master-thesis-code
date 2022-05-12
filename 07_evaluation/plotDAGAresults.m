function [] = plotDAGAresults(evalMixed, evalTarget, evalInterf, ...
    mixedSignal, targetSignal, interfSignal, ...
    dataMixed, dataTarget, dataInterf, AlgorithmParameters, timeVec)

%% subplot DOA histograms

subplot = @(m,n,p) subtightplot(m,n,p,[0.03 0.01], [0.08 0.05], [0.15 0.08]);

figure('position',[2000 50 560 1260]);
s(1) = subplot(3,1,1);
doahistogram(evalTarget, '#D95319')
annotation('textbox','String',"a)",'FontSize',14,'LineStyle','none',...
    'Position',s(1).Position,'FitBoxToText','on')
s(2) = subplot(3,1,2);
doahistogram(evalInterf, '#D95319')
annotation('textbox','String',"b)",'FontSize',14,'LineStyle','none',...
    'Position',s(2).Position,'FitBoxToText','on')
s(3) = subplot(3,1,3);
doahistogram(evalMixed, '#D95319')
annotation('textbox','String',"c)",'FontSize',14,'LineStyle','none',...
    'Position',s(3).Position,'FitBoxToText','on')

%% one plot DOA histograms

figure;
hold on;
doahistogram(evalTarget,'#D95319')
doahistogram(evalInterf,'#EDB120')
doahistogram(evalMixed,'#0072BD')
legend('Target speech signal', 'Interferer speech signal', ...
    'Mixed speech signal')

%% Sectrograms with period glimpse overlay

subplot = @(m,n,p) subtightplot(m,n,p,[0.03 0.01], [0.08 0.05], [0.15 0.08]);
% clear subplot

figure('position',[2000 50 560 860]);
s(1) = subplot(3,1,1);
glimpsespectrogram(targetSignal, timeVec,...
    AlgorithmParameters.Gammatone.samplingRateHz,...
    AlgorithmParameters.Gammatone.nBands,...
    evalTarget, dataTarget.targetSampleIndices, dataTarget.interfSampleIndices)
annotation('textbox','String',"a)",'FontSize',14,'LineStyle','none',...
    'Position',s(1).Position+[0 0.034 0 0],'FitBoxToText','on')
xlabel('Time(s)','Color','none')
set(gca,'xticklabels',[])
s(2) = subplot(3,1,2);
glimpsespectrogram(interfSignal, timeVec,...
    AlgorithmParameters.Gammatone.samplingRateHz,...
    AlgorithmParameters.Gammatone.nBands,...
    evalInterf, dataInterf.targetSampleIndices, dataInterf.interfSampleIndices)
annotation('textbox','String',"b)",'FontSize',14,'LineStyle','none',...
    'Position',s(2).Position+[0 0.034 0 0],'FitBoxToText','on')
xlabel('Time(s)','Color','none')
set(gca,'xticklabels',[])
s(3) = subplot(3,1,3);
glimpsespectrogram(mixedSignal, timeVec,...
    AlgorithmParameters.Gammatone.samplingRateHz,...
    AlgorithmParameters.Gammatone.nBands,...
    evalMixed, dataMixed.targetSampleIndices, dataMixed.interfSampleIndices)
annotation('textbox','String',"c)",'FontSize',14,'LineStyle','none',...
    'Position',s(3).Position+[0 0.034 0 0],'FitBoxToText','on')
set(gca,'xticklabels',[])

xticks(s(3),[0.5 1.5 2.5])
xticklabels(s(3),{'0.5','1.5','2.5'})

end


function [] = doahistogram(evaluation, colorCode)
    histogram(nonzeros(evaluation.azDeg(abs(evaluation.azDeg)<90)),...
        'Normalization', 'probability', 'NumBins', 90, 'FaceColor', colorCode)
%     title(['DOA estimation histogram'])%,['total no. of subband samples: ',num2str(length(testSignal)*nBands)])
%     ylabel({'No. of occurences';'in subband samples'})
    ylabel('Probability of occurrence')
    xlabel('$\hat{\varphi} (^{\circ})$','Interpreter','latex')
    xticks([-90 -60 -30 0 30 60 90])
%     ylim([0 0.25])
    grid on
    set(gca,'FontSize',14);
end

function [] = glimpsespectrogram(testSignal, timeVec, samplingRateHz, ...
    nBands, evaluation, targetSampleIndices, interfSampleIndices)
    meanTestSignal = (testSignal(:,1) + testSignal(:,2))./2;
    winlen = 2000;
    overlap = round(winlen-500);%kaiser5 *0.705);%blackmanharris *0.661);

    spectrogram(meanTestSignal,blackmanharris(winlen),overlap,[],...
        samplingRateHz,'yaxis','power');
    colormap('bone');
    ylim([0, 0.5]);
    caxis([-80 -10])
    hold on;
    for iBand = 1:nBands%7
        scatter(timeVec(targetSampleIndices.L{iBand}), ...
            1e-3.*evaluation.p0SearchRangeFreqVector((evaluation.targetp0DetectedIndexVectors.L{iBand})), ...
            10, 'green','filled');
        scatter(timeVec(interfSampleIndices.L{iBand}), ...
            1e-3.*evaluation.p0SearchRangeFreqVector((evaluation.interfp0DetectedIndexVectors.L{iBand})), ...
            10, 'red','filled');

        scatter(timeVec(targetSampleIndices.R{iBand}), ...
            1e-3.*evaluation.p0SearchRangeFreqVector((evaluation.targetp0DetectedIndexVectors.R{iBand})), ...
            10, 'green','filled');
        scatter(timeVec(interfSampleIndices.R{iBand}), ...
            1e-3.*evaluation.p0SearchRangeFreqVector((evaluation.interfp0DetectedIndexVectors.R{iBand})), ...
            10, 'red','filled');
        legend('target glimpses','interferer glimpses')
    end
    set(gca,'FontSize',14);
end
