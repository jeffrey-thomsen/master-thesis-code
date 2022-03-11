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
AlgorithmParameters.snrCondition    = true;
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
[targetSignal,fsTarget] = audioread('pathological_16.wav');%p298_097.wav');%'intelligence_16.wav');%%'sp01.wav');
[interfSignal,fsInterf] = audioread('peaches_16.wav');%p313_256.wav');%'storylines_16.wav');%%'sp30.wav');

% adjust sampling rate to HRTF
targetSignal = resample(targetSignal,hrtf.Data.SamplingRate,fsTarget);
interfSignal = resample(interfSignal,hrtf.Data.SamplingRate,fsInterf);

% equalize lengths
if length(targetSignal)>length(interfSignal)
    targetSignal = targetSignal(1:length(interfSignal));
elseif length(targetSignal)<length(interfSignal)
    interfSignal = interfSignal(1:length(targetSignal));
end

% equalize levels - SNR 0dB
interfSignal = interfSignal*(sum(abs(targetSignal).^2)./sum(abs(interfSignal).^2));

% convolve with HRTF at different incidence angles
targetSignal = SOFAspat(targetSignal, hrtf, 0, 0);
interfSignal = SOFAspat(interfSignal, hrtf, 60, 0);

% adjust sampling rate to gammatone filterbank
targetSignal = resample(targetSignal,AlgorithmParameters.Gammatone.samplingRateHz,hrtf.Data.SamplingRate);
interfSignal = resample(interfSignal,AlgorithmParameters.Gammatone.samplingRateHz,hrtf.Data.SamplingRate);

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
enhancedSignalMixed = enhancedSignalMixed./max(max(abs(enhancedSignalMixed)));
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


%% Compare coherent periodic samples detected

% number of target samples over all subbands and both channels in each the
% mixed, target and interferer signal
noTSIM = numel(cat(1,dataMixed.targetSampleIndices.L{:},dataMixed.targetSampleIndices.R{:}));
noTSIT = numel(cat(1,dataTarget.targetSampleIndices.L{:},dataTarget.targetSampleIndices.R{:}));
noTSII = numel(cat(1,dataInterf.targetSampleIndices.L{:},dataInterf.targetSampleIndices.R{:}));
% number of interferer samples over all subbands and both channels in each 
% the mixed, target and interferer signal
noISIM = numel(cat(1,dataMixed.interfSampleIndices.L{:},dataMixed.interfSampleIndices.R{:}));
noISIT = numel(cat(1,dataTarget.interfSampleIndices.L{:},dataTarget.interfSampleIndices.R{:}));
noISII = numel(cat(1,dataInterf.interfSampleIndices.L{:},dataInterf.interfSampleIndices.R{:}));

for iBand = 1:30
    % at which sample were target samples identified and with which p0 index?
    eTT = zeros(length(mixedSignal),1);
    eMT = zeros(length(mixedSignal),1);
    eTT(dataTarget.targetSampleIndices.L{iBand}) = evalTarget.targetp0DetectedIndexVectors.L{iBand};
    eMT(dataMixed.targetSampleIndices.L{iBand}) = evalMixed.targetp0DetectedIndexVectors.L{iBand};
    % logical of sample indices where target and mixed signal both identified a target period
    eTT&eMT;
    % mean offset in samples of detected p0
    meanOffset(iBand) = mean(eTT(eTT&eMT) - eMT(eTT&eMT));
    % logicals of samples that were detected either in the mixed OR the target
    eTT(eTT&~eMT);
    eMT(eMT&~eTT);
    % mean detected p0 index
    meanTT(iBand) = mean(p0SearchRangeHz(evalTarget.targetp0DetectedIndexVectors.L{iBand}));
    meanMT(iBand) = mean(p0SearchRangeHz(evalMixed.targetp0DetectedIndexVectors.L{iBand}));
end

figure;
scatter(1:30, meanOffset, 'filled')
title('mean offset of detcted p0 for samples identified as target')
xlabel('subband number')
ylabel('mean offset (samples)')

figure;
hold on;
scatter(1:30, meanTT, 24, 'filled')
scatter(1:30, meanMT, 24, 'filled')
title('mean detcted p0 for samples identified as target')
xlabel('subband number')
ylabel('mean detected p0 (Hz)')
legend('target signal', ' mixed signal')
%% Play mixed signals
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
sound(mixedSignal, AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
sound(enhancedSignalMixed, AlgorithmParameters.Gammatone.samplingRateHz);

%% Play target signals
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
sound(targetSignal, AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
sound(enhancedSignalTarget, AlgorithmParameters.Gammatone.samplingRateHz);

%% Play interferer signals
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
sound(interfSignal, AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
sound(enhancedSignalInterf, AlgorithmParameters.Gammatone.samplingRateHz);

%% Display and save data
% audiowrite(['testInput','.wav'],mixedSignal,AlgorithmParameters.Gammatone.samplingRateHz);
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

    % read position and chosen p0 and azimuth for all detected target and
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

    % proportion of detected periodic samples and samples that passed the
    % coherence mask in azimuth estimation to the toal number of signal
    % samples
    evaluation.periodicityRatio = (numel(nonzeros(cat(1, ...
        p0DetectedIndexVectors.L{:}, p0DetectedIndexVectors.R{:})))) / ...
        (nBands*numel(testSignal));
    evaluation.doaRatio = (numel(nonzeros(evaluation.azDeg))) / ...
        (nBands*length(testSignal));

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
    title('p0 detection histogram')
    xlabel('Period (s)')
    ylabel('No. of occurences in subband samples')

    %% signal spectrogram
    meanTestSignal = (testSignal(:,1) + testSignal(:,2))./2;
    figure;
    winlen = 2000;
    overlap = round(winlen-500);%kaiser5 *0.705);%blackmanharris *0.661);
    spectrogram(meanTestSignal,blackmanharris(winlen),overlap,[],...
        AlgorithmParameters.Gammatone.samplingRateHz,'yaxis');
    ylim([0, 0.5]);
    shading interp;
    title(['signal spectrogram, winlen=',num2str(winlen),', overlap=',num2str(overlap)])
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
            10, 'red','filled');
        scatter(timeVec(targetSampleIndices.R{iBand}), ...
            1e-3.*evaluation.p0SearchRangeFreqVector((evaluation.targetp0DetectedIndexVectors.R{iBand})), ...
            10, 'red','filled');
        scatter(timeVec(interfSampleIndices.L{iBand}), ...
            1e-3.*evaluation.p0SearchRangeFreqVector((evaluation.interfp0DetectedIndexVectors.L{iBand})), ...
            10, 'black','filled');
        scatter(timeVec(interfSampleIndices.R{iBand}), ...
            1e-3.*evaluation.p0SearchRangeFreqVector((evaluation.interfp0DetectedIndexVectors.R{iBand})), ...
            10, 'black','filled');
    end
    
    %% plot azimuth estimation histogram
    figure;
    histogram(nonzeros(evaluation.azDeg(abs(evaluation.azDeg)<90)))
    title('DOA estimation histogram')
    xlabel('azimuth (Degrees)')
    ylabel('No. of occurences in subband samples')

    %% plot position and p0 of detected target and interferer samples
    figure;
    for iBand = 1:30%[3,4,5,10,14,20,25,30]
        subplot(1,2,1)
        title('Left channel')
        xlabel('Time (sec)')
        ylabel('detected frequency (Hz)')
        hold on;
        scatter(timeVec(targetSampleIndices.L{iBand}), ...
            evaluation.p0SearchRangeFreqVector(evaluation.targetp0DetectedIndexVectors.L{iBand}),'red');
        scatter(timeVec(interfSampleIndices.L{iBand}), ...
            evaluation.p0SearchRangeFreqVector(evaluation.interfp0DetectedIndexVectors.L{iBand}),'blue');
        legend('target','interferer')
        ylim(AlgorithmParameters.p0SearchRangeHz)
    
        subplot(1,2,2)
        title('Right channel')
        xlabel('Time (sec)')
        ylabel('detected frequency (Hz)')
        hold on;
        scatter(timeVec(targetSampleIndices.R{iBand}), ...
            evaluation.p0SearchRangeFreqVector(evaluation.targetp0DetectedIndexVectors.R{iBand}),'red');
        scatter(timeVec(interfSampleIndices.R{iBand}), ...
            evaluation.p0SearchRangeFreqVector(evaluation.interfp0DetectedIndexVectors.R{iBand}),'blue');
        legend('target','interferer')
        ylim(AlgorithmParameters.p0SearchRangeHz)
        
        sgtitle(['Processed periodic samples - Gammatone band no. ',num2str(iBand),' fc=',num2str(fc(iBand),'%.0f'),'Hz'])
    end
end