% compute precision and recall rates for DAGA paper
function [precision, recall] = glimpseDistanceMetric(AlgorithmParameters, ...
    dataMixed, dataTarget, dataInterf, centerFreqsHz, mixedSignal)


nBands = AlgorithmParameters.Gammatone.nBands;

metricTargetL = plotDistanceMetric(dataMixed.targetSampleIndices.L,dataTarget.targetSampleIndices.L,...
dataMixed.p0DetectedIndexVectors.L,dataTarget.p0DetectedIndexVectors.L,...
dataMixed.p0SearchRangeSamplesVector, nBands, centerFreqsHz, 'target', length(mixedSignal), false);
metricTargetR = plotDistanceMetric(dataMixed.targetSampleIndices.R,dataTarget.targetSampleIndices.R,...
dataMixed.p0DetectedIndexVectors.R,dataTarget.p0DetectedIndexVectors.R,...
dataMixed.p0SearchRangeSamplesVector, nBands, centerFreqsHz, 'target', length(mixedSignal), false);

metricInterfL = plotDistanceMetric(dataMixed.interfSampleIndices.L,dataInterf.interfSampleIndices.L,...
dataMixed.p0DetectedIndexVectors.L,dataInterf.p0DetectedIndexVectors.L,...
dataMixed.p0SearchRangeSamplesVector, nBands, centerFreqsHz, 'interferer', length(mixedSignal), false);
metricInterfR = plotDistanceMetric(dataMixed.interfSampleIndices.R,dataInterf.interfSampleIndices.R,...
dataMixed.p0DetectedIndexVectors.R,dataInterf.p0DetectedIndexVectors.R,...
dataMixed.p0SearchRangeSamplesVector, nBands, centerFreqsHz, 'interferer', length(mixedSignal), false);


nTotalTruePositivesT = sum(metricTargetL.nTruePositives)+...
                      sum(metricTargetR.nTruePositives);%+...
nTotalTruePositivesI = sum(metricInterfL.nTruePositives)+...
                      sum(metricInterfR.nTruePositives);

nTotalFalsePositivesT = sum(metricTargetL.nFalsePositives)+...
                       sum(metricTargetR.nFalsePositives);%+...
nTotalFalsePositivesI = sum(metricInterfL.nFalsePositives)+...
                       sum(metricInterfR.nFalsePositives);

nTotalFalseNegativesT = sum(metricTargetL.nFalseNegatives)+...
                       sum(metricTargetR.nFalseNegatives);%+...
nTotalFalseNegativesI = sum(metricInterfL.nFalseNegatives)+...
                       sum(metricInterfR.nFalseNegatives);

precision.Target = nTotalTruePositivesT / ...
    (nTotalTruePositivesT + nTotalFalsePositivesT);
precision.Interf = nTotalTruePositivesI / ...
    (nTotalTruePositivesI + nTotalFalsePositivesI);
recall.Target = nTotalTruePositivesT / ...
    (nTotalTruePositivesT + nTotalFalseNegativesT);
recall.Interf = nTotalTruePositivesI / ...
    (nTotalTruePositivesI + nTotalFalseNegativesI);

precision.Total = (nTotalTruePositivesT+nTotalTruePositivesI) / ...
    (nTotalTruePositivesT + nTotalFalsePositivesT + nTotalTruePositivesI + nTotalFalsePositivesI);

end


function Metric = plotDistanceMetric(dataMixedSeparateSampleIndices, ...
    dataSeparateSeparateSampleIndices, dataMixedp0DetectedIndexVectors, ...
    dataSeparatep0DetectedIndexVectors, p0SearchRangeSamplesVector, nBands, ...
    centerFreqsHz, separateString, lengthOfSignal, Plotting)


truePositiveLogical = cell(1,nBands);
truePositiveSampleIndices = cell(1,nBands);

falsePositiveLogical = cell(1,nBands);
falsePositiveSampleIndices = cell(1,nBands);

falseNegativeLogical = cell(1,nBands);

Metric.nTruePositives  = zeros(1,nBands);
Metric.nFalsePositives = zeros(1,nBands);
Metric.nFalseNegatives = zeros(1,nBands);

Metric.precision = zeros(1,nBands);
Metric.recall = zeros(1,nBands);

mixedTargetP0Indices = cell(1,nBands);
mixedTargetP0ValuesSamples = cell(1,nBands);
targetTargetP0Indices = cell(1,nBands);
targetTargetP0ValuesSamples = cell(1,nBands);
Metric.nSharedGlimpses = zeros(1,nBands);
Metric.nSharedPeriods = zeros(1,nBands);
distanceOfSharedGlimpsesInSamples = cell(1,nBands);
meanAbsDistanceOfTargetPeriodsInSamples = zeros(1,nBands);
Metric.nOverallGlimpses = zeros(1,nBands);

for iBand = 1:nBands
if ~isempty(dataMixedSeparateSampleIndices{iBand})
% find which target glimpses in the mixed signal are also in the target
% signal
truePositiveLogical{iBand} = ismember(dataMixedSeparateSampleIndices{iBand},dataSeparateSeparateSampleIndices{iBand});
% count them
Metric.nTruePositives(iBand) = nnz(truePositiveLogical{iBand});

% get the sample indices of those glimpses
truePositiveSampleIndices{iBand} = dataMixedSeparateSampleIndices{iBand}(truePositiveLogical{iBand});
% get the detected periods for those glimpses in the mixed signal
mixedTargetP0Indices{iBand} = dataMixedp0DetectedIndexVectors{iBand}(truePositiveSampleIndices{iBand});
mixedTargetP0ValuesSamples{iBand} = p0SearchRangeSamplesVector(mixedTargetP0Indices{iBand});
% get the detected periods for those glimpses in the target signal
targetTargetP0Indices{iBand} = dataSeparatep0DetectedIndexVectors{iBand}(truePositiveSampleIndices{iBand});
targetTargetP0ValuesSamples{iBand} = p0SearchRangeSamplesVector(targetTargetP0Indices{iBand});

% find which target glimpses in the mixed signal are not in the target signal
falsePositiveLogical{iBand} = ~truePositiveLogical{iBand};
% count them
Metric.nFalsePositives(iBand) = nnz(falsePositiveLogical{iBand});
% get the sample indices of those glimpses
falsePositiveSampleIndices{iBand} = dataMixedSeparateSampleIndices{iBand}(falsePositiveLogical{iBand});

% find which target glimpses in the target signal are not in the mixed signal
falseNegativeLogical{iBand} = ~ismember(dataSeparateSeparateSampleIndices{iBand}, dataMixedSeparateSampleIndices{iBand});
% count them
Metric.nFalseNegatives(iBand) = nnz(falseNegativeLogical{iBand});

Metric.precision(iBand) = Metric.nTruePositives(iBand) / ...
    (Metric.nTruePositives(iBand)+Metric.nFalsePositives(iBand));
Metric.recall(iBand) = Metric.nTruePositives(iBand) / ...
    (Metric.nTruePositives(iBand)+Metric.nFalseNegatives(iBand));



Metric.nSharedGlimpses(iBand) = numel(mixedTargetP0ValuesSamples{iBand});
Metric.nSharedPeriods(iBand) = nnz(mixedTargetP0ValuesSamples{iBand}==targetTargetP0ValuesSamples{iBand});
distanceOfSharedGlimpsesInSamples{iBand} = mixedTargetP0ValuesSamples{iBand}-targetTargetP0ValuesSamples{iBand};
meanAbsDistanceOfTargetPeriodsInSamples(iBand) = mean(abs(distanceOfSharedGlimpsesInSamples{iBand}));

Metric.nOverallGlimpses(iBand) = numel(dataMixedSeparateSampleIndices{iBand});
% if iBand==5
% figure;
% hold on
% plot(mixedTargetP0ValuesSamples{iBand})
% plot(targetTargetP0ValuesSamples{iBand})
% title(['Shared ',separateString,' glimpses - Gammatone band no. ',num2str(iBand)],['fc=',num2str(centerFreqsHz(iBand),'%.0f'),'Hz'])
% legend('Mixed signal', [separateString,' signal'])
% xlabel('x')
% ylabel('y')
% end
end
end
if Plotting
    % figure;
    % scatter(centerFreqsHz,nSharedPeriods./nSharedGlimpses)
    % xlabel('gammatone center frequency (Hz)')
    % ylabel('ratio (1)')
    % title({['ratio of no. of shared ',separateString,' glimpses with identical period'],['to overall no. of shared ',separateString,' glimpses']})

    figure;
    hold on;
    bar(centerFreqsHz, Metric.nOverallGlimpses')
    bar(centerFreqsHz, Metric.nSharedGlimpses')
    bar(centerFreqsHz, Metric.nSharedPeriods')
    % scatter(centerFreqsHz,nSharedTargetPeriods)
    % scatter(centerFreqsHz,nSharedTargetSamples)
    xlabel('gammatone center frequency (Hz)')
    ylabel('no. of glimpses')
    legend(['overall no. of ',separateString,' glimpses in mixed signal'], ['no. of shared ',separateString,' glimpses'], ['no. of shared ',separateString,' glimpses with same p_0'])
    xlim([0 1500])
    title(['Comparing ',separateString,' glimpses in mixed and ',separateString,' signals'],['overall no. of samples per subband: ',num2str(lengthOfSignal)])
    set(gca,'FontSize',14);
    
    figure;
    hold on;
    % bar(centerFreqsHz,[Metric.precision; Metric.recall],1.6)
    bar(centerFreqsHz,Metric.precision)
    % bar(centerFreqsHz,Metric.recall)
    % scatter(centerFreqsHz, Metric.precision,25, 'red','filled')
    % scatter(centerFreqsHz, Metric.recall, 25, 'blue','filled')
    xlabel('gammatone center frequency (Hz)')
    ylabel('ratio (1)')
    legend([separateString,' precision'])%, [separateString,' recall'])
    xlim([0 1500])
    ylim([0 1])
    title(['Comparing ',separateString,' glimpses in mixed and ',separateString,' signals'])
    set(gca,'FontSize',14);
    
    % figure;
    % hold on
    % bar(centerFreqsHz,meanAbsDistanceOfTargetPeriodsInSamples)
    % % scatter(centerFreqsHz,nSharedTargetPeriods)
    % % scatter(centerFreqsHz,nSharedTargetSamples)
    % xlabel('gammatone center frequency (Hz)')
    % ylabel('mean distance (samples)')
    % title({['Mean absolute distance of detected periods'],['of shared ',separateString,' glimpses']})
    % xlim([0 1000])
end

end
