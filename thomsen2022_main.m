% This is the MATLAB code for the master thesis "Speech enhacment for
% hearing aids using a combination of binaural and periodicity features" by
% Jeffrey Thomsen, written at the CvO University of Oldenburg, Germany, in
% 2022.
% This is the main script, the running of which will allow you to define
% custom simulation, algorithm and evaluation parameters, or rerun a 
% simulation from a log file. The simulation is executed and the data and
% evaluation results are displayed and stored accordingly.

%% Start/Simulation interface
clear
% msgbox('Hi, you are running the Thomsen2022 speech enhancement algorithm!')

%% Preparation

% initialize structs
% TestSignalParameters = struct;
% TargetAngleParameters = struct;
% BlockFeedingParameters = struct;
AlgorithmParameters = AlgorithmParametersConstructor();

% Alterations of simulation parameters
AlgorithmParameters.p0SearchRangeHz = [100 350];
AlgorithmParameters.ChenP0Detection = false;
AlgorithmParameters.coherenceMask   = true;
AlgorithmParameters.azimuthPooling  = false;
AlgorithmParameters.snrCondition    = false;
AlgorithmParameters.DOAProcessing   = true;
AlgorithmParameters.Cancellation    = true;
AlgorithmParameters.Enhancement     = true;
AlgorithmParameters.snrLPFilterTau = 0.04;
AlgorithmParameters.ivsThreshold = 0.98; % Dietz2011
AlgorithmParameters.nCyclesTau = 5; % Dietz2011

AlgorithmParameters.snrThresholdInDb = 10;

% load HRTF
hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
AlgorithmParameters.Gammatone.samplingRateHz = hrtf.Data.SamplingRate;

% generate/load IPD-to-azimuth mapping function
tic
load('2022-03-04_itd_lookuptable_annotated.mat');
lookuptable = lookuptable.lookuptable;
AlgorithmParameters.lookuptable = lookuptable;
% lookuptable = interauralToAzimuthLookuptable(hrtf, AlgorithmParameters,...
%     'intelligence_16.wav', 'storylines_16.wav', 'shellshock_16.wav', ...
%     'peaches_16.wav', 'p360_253.wav', 'p313_256.wav', 'p298_097.wav', ...
%     'p237_079.wav', 'p234_003.wav', '2078-142845-0002.flac');
% interauralToAzimuthLookuptable(hrtf, AlgorithmParameters,'intelligence_16.wav', 'storylines_16.wav', 'shellshock_16.wav', 'peaches_16.wav', 'housewives_16.wav', 'necessity_16.wav', 'prison_16.wav', 'butterscotch_16.wav', 'bigtips_16.wav', 'pathological_16.wav');
toc


% set sampling frequency to at least 2x upper gammatone centre frequency
AlgorithmParameters.Gammatone.samplingRateHz = ...
    ceil(2*AlgorithmParameters.Gammatone.fHigh);

% generate gammatone filterbank
[AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);

centerFreqsHz = AlgorithmStates.L.GammatoneStates.analyzer.center_frequencies_hz;

%% Generate test signal

% load clean speech signals
[targetSignal, fsTarget] = audioread('pathological_16.wav');%p298_097.wav');%'intelligence_16.wav');%%'sp01.wav');
[interfSignal, fsInterf] = audioread('peaches_16.wav');%p313_256.wav');%'storylines_16.wav');%%'sp30.wav');

% adjust sampling rate to HRTF
targetSignal = resample(targetSignal, hrtf.Data.SamplingRate, fsTarget);
interfSignal = resample(interfSignal, hrtf.Data.SamplingRate, fsInterf);

% equalize lengths
if length(targetSignal) > length(interfSignal)
    targetSignal = targetSignal(1:length(interfSignal));
elseif length(targetSignal) < length(interfSignal)
    interfSignal = interfSignal(1:length(targetSignal));
end

% equalize levels - SNR 0dB
% interfSignal = interfSignal*(sum(abs(targetSignal).^2)./sum(abs(interfSignal).^2));

% convolve with HRTF at different incidence angles
targetSignal = SOFAspat(targetSignal, hrtf, 0, 0);
interfSignal = SOFAspat(interfSignal, hrtf, 60, 0);

% adjust sampling rate to gammatone filterbank
targetSignal = resample(targetSignal, ...
    AlgorithmParameters.Gammatone.samplingRateHz, hrtf.Data.SamplingRate);
interfSignal = resample(interfSignal, ...
    AlgorithmParameters.Gammatone.samplingRateHz, hrtf.Data.SamplingRate);

% equalize levels - SIR 0dB
targetSignal = targetSignal / std(targetSignal(:));
interfSignal = interfSignal / std(interfSignal(:));

% add and normalize test signal
mixedSignal = targetSignal + interfSignal;
mixedSignal = mixedSignal./max(max(mixedSignal));

% testSignal = testSignalGenerator;
mixedSignal = testInputSignal(mixedSignal);

% time vector for plotting
dt = 1/AlgorithmParameters.Gammatone.samplingRateHz;
timeVec = dt:dt:dt*length(mixedSignal);

%% Process mixed speech
% tic
% processedSignal = blockFeedingRoutine(testSignal,BlockFeedingParameters,...
%     AlgorithmParameters, AlgorithmStates);
% toc
tic
[enhancedSignalMixed, AlgorithmStatesMixed, dataMixed] = ...
    speechEnhancement(mixedSignal, AlgorithmParameters, AlgorithmStates);
toc
enhancedSignalMixed = enhancedSignalMixed ./ ...
    max(max(abs(enhancedSignalMixed)));
%% Evaluate mixed speech

evalMixed = evaluateAlgorithm(dataMixed.p0DetectedIndexVectors, AlgorithmParameters, ...
    dataMixed.p0SearchRangeSamplesVector, mixedSignal, timeVec, dataMixed.ivsMaskCells, ...
    dataMixed.azimuthDegCells, dataMixed.targetSampleIndices, dataMixed.interfSampleIndices, centerFreqsHz);

%% Process target speech
tic
[enhancedSignalTarget, AlgorithmStates, dataTarget] = ...
    speechEnhancement(targetSignal, AlgorithmParameters, AlgorithmStates);
toc
enhancedSignalTarget = enhancedSignalTarget./max(max(abs(enhancedSignalTarget)));
%% Evaluate target speech

evalTarget = evaluateAlgorithm(dataTarget.p0DetectedIndexVectors, AlgorithmParameters, ...
    dataTarget.p0SearchRangeSamplesVector, targetSignal, timeVec, dataTarget.ivsMaskCells, ...
    dataTarget.azimuthDegCells, dataTarget.targetSampleIndices, dataTarget.interfSampleIndices, centerFreqsHz);

p0SearchRangeHz = evalMixed.p0SearchRangeFreqVector;
%% Process interferer speech
tic
[enhancedSignalInterf, AlgorithmStatesInterf, dataInterf] = ...
    speechEnhancement(interfSignal, AlgorithmParameters, AlgorithmStates);
toc
enhancedSignalInterf = enhancedSignalInterf./max(max(abs(enhancedSignalInterf)));
%% Evaluate interferer speech

evalInterf = evaluateAlgorithm(dataInterf.p0DetectedIndexVectors, AlgorithmParameters, ...
    dataInterf.p0SearchRangeSamplesVector, interfSignal, timeVec, dataInterf.ivsMaskCells, ...
    dataInterf.azimuthDegCells, dataInterf.targetSampleIndices, dataInterf.interfSampleIndices, centerFreqsHz);



% %% Play mixed signals
% box1 = msgbox('Play original signal (Ensure volume is adequately set)');
% waitfor(box1);
% sound(mixedSignal, AlgorithmParameters.Gammatone.samplingRateHz);
% box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
%     ' section of script (adjusting filter variables if necessary)']);
% waitfor(box2);
% sound(enhancedSignalMixed, AlgorithmParameters.Gammatone.samplingRateHz);
% 
% %% Play target signals
% box1 = msgbox('Play original signal (Ensure volume is adequately set)');
% waitfor(box1);
% sound(targetSignal, AlgorithmParameters.Gammatone.samplingRateHz);
% box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
%     ' section of script (adjusting filter variables if necessary)']);
% waitfor(box2);
% sound(enhancedSignalTarget, AlgorithmParameters.Gammatone.samplingRateHz);
% 
% %% Play interferer signals
% box1 = msgbox('Play original signal (Ensure volume is adequately set)');
% waitfor(box1);
% sound(interfSignal, AlgorithmParameters.Gammatone.samplingRateHz);
% box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
%     ' section of script (adjusting filter variables if necessary)']);
% waitfor(box2);
% sound(enhancedSignalInterf, AlgorithmParameters.Gammatone.samplingRateHz);

%% Save data
audiowrite(['testInput','.wav'],mixedSignal,AlgorithmParameters.Gammatone.samplingRateHz);
% audiowrite(['testJeffrey','.wav'],enhancedMixedSignal,AlgorithmParameters.Gammatone.samplingRateHz);

MetaData = struct;
[~,MetaData.gitCommitHash]=system('git rev-parse --short HEAD');
gitURL = 'https://github.com/jeffrey-thomsen/master-thesis-code/commit/';
MetaData.gitCommitURL = strcat(gitURL, MetaData.gitCommitHash);
MetaData.date = datetime('now','Format','yyyy-MM-dd''_''HH-mm');
dateString = string(MetaData.date);

filename = string(dateString+'_simulation_data.mat');
save(filename,...
'AlgorithmParameters',...
'AlgorithmStates', 'AlgorithmStatesInterf', 'AlgorithmStatesMixed',...
'centerFreqsHz',...
'dataInterf', 'dataMixed', 'dataTarget',...
'enhancedSignalInterf', 'enhancedSignalMixed', 'enhancedSignalTarget',...
'evalInterf', 'evalMixed', 'evalTarget',...
'interfSignal', 'mixedSignal', 'targetSignal',...
'lookuptable',...
'MetaData')

%% definitive glimpse distance metric
nBands = AlgorithmParameters.Gammatone.nBands;

metricTargetL = plotDistanceMetric(dataMixed.targetSampleIndices.L,dataTarget.targetSampleIndices.L,...
dataMixed.p0DetectedIndexVectors.L,dataTarget.p0DetectedIndexVectors.L,...
dataMixed.p0SearchRangeSamplesVector, nBands, centerFreqsHz, 'target', length(mixedSignal));
metricTargetR = plotDistanceMetric(dataMixed.targetSampleIndices.R,dataTarget.targetSampleIndices.R,...
dataMixed.p0DetectedIndexVectors.R,dataTarget.p0DetectedIndexVectors.R,...
dataMixed.p0SearchRangeSamplesVector, nBands, centerFreqsHz, 'target', length(mixedSignal));

metricInterfL = plotDistanceMetric(dataMixed.interfSampleIndices.L,dataInterf.interfSampleIndices.L,...
dataMixed.p0DetectedIndexVectors.L,dataInterf.p0DetectedIndexVectors.L,...
dataMixed.p0SearchRangeSamplesVector, nBands, centerFreqsHz, 'interferer', length(mixedSignal));
metricInterfR = plotDistanceMetric(dataMixed.interfSampleIndices.R,dataInterf.interfSampleIndices.R,...
dataMixed.p0DetectedIndexVectors.R,dataInterf.p0DetectedIndexVectors.R,...
dataMixed.p0SearchRangeSamplesVector, nBands, centerFreqsHz, 'interferer', length(mixedSignal));


nTotalTruePositivesT = sum(metricTargetL.nTruePositives)+...
                      sum(metricTargetR.nTruePositives);%+...
nTotalTruePositivesI =                       sum(metricInterfL.nTruePositives)+...
                      sum(metricInterfR.nTruePositives);

nTotalFalsePositivesT = sum(metricTargetL.nFalsePositives)+...
                       sum(metricTargetR.nFalsePositives);%+...
nTotalFalsePositivesI =                        sum(metricInterfL.nFalsePositives)+...
                       sum(metricInterfR.nFalsePositives);

nTotalFalseNegativesT = sum(metricTargetL.nFalseNegatives)+...
                       sum(metricTargetR.nFalseNegatives);%+...
nTotalFalseNegativesI =                        sum(metricInterfL.nFalseNegatives)+...
                       sum(metricInterfR.nFalseNegatives);

precision = nTotalTruePositives / ...
    (nTotalTruePositives + nTotalFalsePositives);
recall = nTotalTruePositives / ...
    (nTotalTruePositives + nTotalFalseNegatives);
%%
plotGlimpseOverlay(targetSignal, AlgorithmParameters.Gammatone.samplingRateHz, timeVec, nBands,...
    dataTarget.targetSampleIndices, dataMixed.targetSampleIndices, p0SearchRangeHz,...
    evalTarget.targetp0DetectedIndexVectors, evalMixed.targetp0DetectedIndexVectors, 'target')


%% Auxiliary functions

function evaluation = evaluateAlgorithm(p0DetectedIndexVectors, ...
    AlgorithmParameters, p0SearchRangeSamplesVector, testSignal, timeVec, ...
    ivsMask, azimuthDegCells, targetSampleIndices, interfSampleIndices, fc)

    nBands = AlgorithmParameters.Gammatone.nBands;
        
    % map detected sample indices to p0 values
    for iBand = 1:nBands
        evaluation.targetp0DetectedIndexVectors.L{iBand} = ...
            p0DetectedIndexVectors.L{iBand}(targetSampleIndices.L{iBand});
        evaluation.interfp0DetectedIndexVectors.L{iBand} = ...
            p0DetectedIndexVectors.L{iBand}(interfSampleIndices.L{iBand});
        evaluation.targetp0DetectedIndexVectors.R{iBand} = ...
            p0DetectedIndexVectors.R{iBand}(targetSampleIndices.R{iBand});
        evaluation.interfp0DetectedIndexVectors.R{iBand} = ...
            p0DetectedIndexVectors.R{iBand}(interfSampleIndices.R{iBand});
    end

    % read detected azimuth and position for all detected target and
    %  interferer samples
    evaluation.azDeg = zeros(nBands,length(ivsMask{1}));
    for iBand = 1:nBands
        if AlgorithmParameters.coherenceMask
            evaluation.azDeg(iBand,ivsMask{iBand}) = azimuthDegCells{iBand};
        else
            evaluation.azDeg(iBand,:) = azimuthDegCells{iBand};
        end
        evaluation.coherentPeriodicComponentLogicalVector.L(iBand,:) = ...
            ivsMask{iBand} & p0DetectedIndexVectors.L{iBand};
        evaluation.coherentPeriodicComponentLogicalVector.R(iBand,:) = ...
            ivsMask{iBand} & p0DetectedIndexVectors.R{iBand};
    end

    % ratio of detected periodic samples and samples that passed the
    % coherence mask in azimuth estimation and coherent periodic samples
    % to the toal number of subband signal samples
    evaluation.periodicityRatio = (nnz(cat(1, ...
        p0DetectedIndexVectors.L{:}, p0DetectedIndexVectors.R{:}))) / ...
        (nBands*numel(testSignal));
    evaluation.doaRatio = (nnz(evaluation.azDeg)) / ...
        (nBands*length(testSignal));
    evaluation.coherentPeriodicRatio = ...
    (nnz(evaluation.coherentPeriodicComponentLogicalVector.L) + ...
     nnz(evaluation.coherentPeriodicComponentLogicalVector.L)) /...
     numel(evaluation.coherentPeriodicComponentLogicalVector.L);

    %% p0 detection histogram
    [evaluation.GC,evaluation.GR] = groupcounts(cat(1,p0DetectedIndexVectors.L{:},...
        p0DetectedIndexVectors.R{:}));
    evaluation.GR = nonzeros(evaluation.GR);%+p0SearchRangeSamplesVector(1)-1;
    evaluation.GC = evaluation.GC(end-length(evaluation.GR)+1:end);
    evaluation.p0SearchRangeSecondsVector = ...
        p0SearchRangeSamplesVector./AlgorithmParameters.Gammatone.samplingRateHz;
    evaluation.p0SearchRangeFreqVector = ...
        AlgorithmParameters.Gammatone.samplingRateHz./p0SearchRangeSamplesVector;
    
    figure;
    hBar = bar(evaluation.p0SearchRangeSecondsVector(evaluation.GR),...
        evaluation.GC);
%     set(gca,'XScale','log');
    title('p0 detection histogram',['total no. of subband samples: ',num2str(numel(testSignal)*nBands)])
    xlabel('Period (s)')
    ylabel('No. of occurences in subband samples')
    set(gca,'FontSize',14);
    %% signal spectrogram
    meanTestSignal = (testSignal(:,1) + testSignal(:,2))./2;
    winlen = 2000;
    overlap = round(winlen-500);%kaiser5 *0.705);%blackmanharris *0.661);
    
    figure;
    spectrogram(meanTestSignal,blackmanharris(winlen),overlap,[],...
        AlgorithmParameters.Gammatone.samplingRateHz,'yaxis');
    colormap('bone');
    ylim([0, 0.5]);
    % title(['signal spectrogram, winlen=',num2str(winlen),', overlap=',num2str(overlap)])
    title('Spectrogram with glimpse overlay')
    hold on;
    for iBand = 1:nBands%7
%         scatter(timeVec(logical(p0DetectedIndexVectors.L{iBand})), ...
%             1e-3.*evaluation.p0SearchRangeFreqVector(nonzeros(p0DetectedIndexVectors.L{iBand})), ...
%             5, 'red','filled');
%         scatter(timeVec(logical(p0DetectedIndexVectors.R{iBand})), ...
%             1e-3.*evaluation.p0SearchRangeFreqVector(nonzeros(p0DetectedIndexVectors.R{iBand})), ...
%             5, 'red','filled');
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
    
    figure;
    spectrogram(meanTestSignal,blackmanharris(winlen),overlap,[],...
        AlgorithmParameters.Gammatone.samplingRateHz,'yaxis');
    colormap('bone');
    ylim([0, 0.5]);
%     title(['signal spectrogram, winlen=',num2str(winlen),', overlap=',num2str(overlap)])
    title('Spectrogram with glimpse overlay')
    hold on;
    for iBand = 1:nBands%7
        scatter(timeVec(logical(p0DetectedIndexVectors.L{iBand})), ...
            1e-3.*evaluation.p0SearchRangeFreqVector(nonzeros(p0DetectedIndexVectors.L{iBand})), ...
            10, 'green','filled');
        scatter(timeVec(logical(p0DetectedIndexVectors.R{iBand})), ...
            1e-3.*evaluation.p0SearchRangeFreqVector(nonzeros(p0DetectedIndexVectors.R{iBand})), ...
            10, 'green','filled');
        legend('periodic glimpses')
    end
    set(gca,'FontSize',14);

    %% plot azimuth estimation histogram
    figure;
    histogram(nonzeros(evaluation.azDeg(abs(evaluation.azDeg)<90)))
    title(['DOA estimation histogram'])%,['total no. of subband samples: ',num2str(length(testSignal)*nBands)])
    ylabel('No. of occurences in subband samples')
    xlabel('azimuth (Degrees)')
    xticks([-90 -60 -30 0 30 60 90])
    set(gca,'FontSize',14);


    %% plot position and p0 of detected target and interferer samples
%     iStart = 1;
%     iEnd = 30;
%     figure;
%     for iBand = iStart:iEnd%[3,4,5,10,14,20,25,30]
%         subplot(1,2,1)
%         title('Left channel')
%         xlabel('Time (s)')
%         ylabel('Detected period (s)')
%         hold on;
%         scatter(timeVec(targetSampleIndices.L{iBand}), ...
%             evaluation.p0SearchRangeSecondsVector(evaluation.targetp0DetectedIndexVectors.L{iBand}),'red');
%         scatter(timeVec(interfSampleIndices.L{iBand}), ...
%             evaluation.p0SearchRangeSecondsVector(evaluation.interfp0DetectedIndexVectors.L{iBand}),'blue');
%         legend('target glimpses','interferer glimpses')
%         ylim(flip(1./AlgorithmParameters.p0SearchRangeHz))
%     
%         subplot(1,2,2)
%         title('Right channel')
%         xlabel('Time (s)')
%         ylabel('Detected period (s)')
%         hold on;
%         scatter(timeVec(targetSampleIndices.R{iBand}), ...
%             evaluation.p0SearchRangeSecondsVector(evaluation.targetp0DetectedIndexVectors.R{iBand}),'red');
%         scatter(timeVec(interfSampleIndices.R{iBand}), ...
%             evaluation.p0SearchRangeSecondsVector(evaluation.interfp0DetectedIndexVectors.R{iBand}),'blue');
%         legend('target glimpses','interferer glimpses')
%         ylim(flip(1./AlgorithmParameters.p0SearchRangeHz))
%     end
%     sgtitle({['Processed glimpses - Gammatone band no. ',num2str(iStart),'-',num2str(iEnd)],['fc=',num2str(fc(iStart),'%.0f'),'-',num2str(fc(iEnd),'%.0f'),'Hz']})
end



function Metric = plotDistanceMetric(dataMixedSeparateSampleIndices, dataSeparateSeparateSampleIndices, ...
    dataMixedp0DetectedIndexVectors, dataSeparatep0DetectedIndexVectors, ...
    p0SearchRangeSamplesVector, nBands, centerFreqsHz, separateString, lengthOfSignal)


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



function [] = plotGlimpseOverlay(testSignal, samplingRateHz, timeVec, nBands,...
    separateSampleIndices, mixedSampleIndices, p0SearchRangeFreqVector,...
    separatep0DetectedIndexVectors, mixedp0DetectedIndexVectors, separateString)
    meanTestSignal = (testSignal(:,1) + testSignal(:,2))./2;
    figure;
    winlen = 2000;
    overlap = round(winlen-500);%kaiser5 *0.705);%blackmanharris *0.661);
%     spectrogram(meanTestSignal,blackmanharris(winlen),overlap,[],...
%         samplingRateHz,'yaxis');
    ylim([0, 0.5]);
    shading interp;
    title(['signal spectrogram, winlen=',num2str(winlen),', overlap=',num2str(overlap)])
    hold on;
    for iBand = 1:nBands%7
        scatter(timeVec(separateSampleIndices.L{iBand}), ...
            1e-3.*p0SearchRangeFreqVector(separatep0DetectedIndexVectors.L{iBand}), ...
            15, 'red','filled');
        scatter(timeVec(mixedSampleIndices.L{iBand}), ...
            1e-3.*p0SearchRangeFreqVector(mixedp0DetectedIndexVectors.L{iBand}), ...
            10, 'black','filled');
        
        scatter(timeVec(separateSampleIndices.R{iBand}), ...
            1e-3.*p0SearchRangeFreqVector(separatep0DetectedIndexVectors.R{iBand}), ...
            15, 'red','filled');
        scatter(timeVec(mixedSampleIndices.R{iBand}), ...
            1e-3.*p0SearchRangeFreqVector(mixedp0DetectedIndexVectors.R{iBand}), ...
            10, 'black','filled');
        legend(['glimpses from ',separateString,' signal'], 'glimpses from mixed signal')
    end
end