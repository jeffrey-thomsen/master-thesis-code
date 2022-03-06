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
TestSignalParameters = struct;
TargetAngleParameters = struct;
BlockFeedingParameters = struct;
AlgorithmParameters = AlgorithmParametersConstructor();

% Alterations
AlgorithmParameters.p0SearchRangeHz = [100 350];
AlgorithmParameters.DOAProcessing = true;

% load HRTF
hrtf = SOFAload('HRIR_KEMAR_DV0001_3.sofa',[5 2],'R');
AlgorithmParameters.Gammatone.samplingRateHz = hrtf.Data.SamplingRate;

% training data strings


% generate/load IPD-to-azimuth mapping function
tic
lookuptable = load('2022-03-04_itd_lookuptable.mat');
lookuptable = lookuptable.lookuptable;
AlgorithmParameters.lookuptable = lookuptable;
%     interauralToAzimuthLookuptable(hrtf, AlgorithmParameters,...
%     'intelligence_16.wav', 'storylines_16.wav');
toc

% set sampling frequency to at least 2x upper gammatone centre frequency
AlgorithmParameters.Gammatone.samplingRateHz = ...
    ceil(2*AlgorithmParameters.Gammatone.fHigh);

% generate gammatone filterbank
[AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);

%% Generate test signal

% load clean speech signals
[targetSignal,fsTarget]=audioread('intelligence_16.wav');%'p298_097.wav');%'sp01.wav');
[interfSignal,fsInterf]=audioread('storylines_16.wav');%'p313_256.wav');%'sp30.wav');

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

%% Block-feeding routine - Speech enhancement algorithm
% tic
% processedSignal = blockFeedingRoutine(testSignal,BlockFeedingParameters,...
%     AlgorithmParameters, AlgorithmStates);
% toc
tic
[enhancedSignalMixed, AlgorithmStatesMixed, dataMixed] = ...
    speechEnhancement(mixedSignal, AlgorithmParameters, AlgorithmStates);
toc
enhancedSignalMixed = enhancedSignalMixed./max(max(abs(enhancedSignalMixed)));
%% Evaluation of mixed speech

evalMixed = evaluateAlgorithm(dataMixed.p0DetectedIndexVectors, AlgorithmParameters, ...
    dataMixed.p0SearchRangeSamplesVector, mixedSignal, dataMixed.ivsMaskCells, ...
    dataMixed.azimuthDegCells, dataMixed.targetSampleIndices, dataMixed.interfSampleIndices);

%% Process target speech
tic
[enhancedSignalTarget, AlgorithmStates, dataTarget] = ...
    speechEnhancement(targetSignal, AlgorithmParameters, AlgorithmStates);
toc
enhancedSignalTarget = enhancedSignalTarget./max(max(abs(enhancedSignalTarget)));
%% Evaluation of target speech

evalTarget = evaluateAlgorithm(dataTarget.p0DetectedIndexVectors, AlgorithmParameters, ...
    dataTarget.p0SearchRangeSamplesVector, mixedSignal, dataTarget.ivsMaskCells, ...
    dataTarget.azimuthDegCells, dataTarget.targetSampleIndices, dataTarget.interfSampleIndices);

%% Process interferer speech
tic
[enhancedSignalInterf, AlgorithmStatesInterf, dataInterf] = ...
    speechEnhancement(interfSignal, AlgorithmParameters, AlgorithmStates);
toc
enhancedSignalInterf = enhancedSignalInterf./max(max(abs(enhancedSignalInterf)));
%% Evaluation of interferer speech

evalInterf = evaluateAlgorithm(dataInterf.p0DetectedIndexVectors, AlgorithmParameters, ...
    dataInterf.p0SearchRangeSamplesVector, mixedSignal, dataInterf.ivsMaskCells, ...
    dataInterf.azimuthDegCells, dataInterf.targetSampleIndices, dataInterf.interfSampleIndices);


%% Compare coherent periodic samples detected

%
noTSIM = numel(cat(1,dataMixed.targetSampleIndices.L{:},dataMixed.targetSampleIndices.R{:}));
noTSIT = numel(cat(1,dataTarget.targetSampleIndices.L{:},dataTarget.targetSampleIndices.R{:}));
noTSII = numel(cat(1,dataInterf.targetSampleIndices.L{:},dataInterf.targetSampleIndices.R{:}));

noISIM = numel(cat(1,dataMixed.interfSampleIndices.L{:},dataMixed.interfSampleIndices.R{:}));
noISIT = numel(cat(1,dataTarget.interfSampleIndices.L{:},dataTarget.interfSampleIndices.R{:}));
noISII = numel(cat(1,dataInterf.interfSampleIndices.L{:},dataInterf.interfSampleIndices.R{:}));



%% Play signals
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
sound(mixedSignal, AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
sound(enhancedSignalMixed, AlgorithmParameters.Gammatone.samplingRateHz);

%% Play signals target
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
sound(targetSignal, AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
sound(enhancedSignalTarget, AlgorithmParameters.Gammatone.samplingRateHz);

%% Play signals interferer
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

%% Auxiliary functions
function evaluation = evaluateAlgorithm(p0DetectedIndexVectors, ...
    AlgorithmParameters, p0SearchRangeSamplesVector, testSignal, ...
    ivsMask, azimuthDegCells, targetSampleIndices, interfSampleIndices)

    [evaluation.GC,evaluation.GR] = groupcounts(cat(1,p0DetectedIndexVectors.L{:},...
        p0DetectedIndexVectors.R{:}));
    evaluation.GR = nonzeros(evaluation.GR);%+p0SearchRangeSamplesVector(1)-1;
    evaluation.GC = evaluation.GC(end-length(evaluation.GR)+1:end);
    evaluation.p0SearchRangeFreqVector = ...
        AlgorithmParameters.Gammatone.samplingRateHz./p0SearchRangeSamplesVector;
    
%     figure;
%     title('p0 detection histogram')
%     xlabel('Frequency (Hz)')
%     ylabel('No. of occurences in subband samples')
%     hBar = bar(evaluation.p0SearchRangeFreqVector(evaluation.GR),...
%         evaluation.GC);
%     set(gca,'XScale','log');

    meanTestSignal = (testSignal(:,1) + testSignal(:,2))./2;
%     figure;
%     title('signal spectrogram')
%     spectrogram(meanTestSignal,hamming(1000),[],[],...
%         AlgorithmParameters.Gammatone.samplingRateHz);xlim([0, 0.5]);
    
    evaluation.azDeg = zeros(30,length(ivsMask{1}));
    for iBand = 1:30
        evaluation.azDeg(iBand,ivsMask{iBand}) = azimuthDegCells{iBand};
        evaluation.coherentPeriodicComponentLogicalVector.L(iBand,:) = ...
            ivsMask{iBand} & p0DetectedIndexVectors.L{iBand};
        evaluation.coherentPeriodicComponentLogicalVector.R(iBand,:) = ...
            ivsMask{iBand} & p0DetectedIndexVectors.R{iBand};
    end
    
%     figure;
%     title('DOA estimation histogram')
%     xlabel('azimuth (Degrees)')
%     ylabel('No. of occurences in subband samples')
%     histogram(nonzeros(evaluation.azDeg(abs(evaluation.azDeg)<90)))
    
    evaluation.periodicityRatio = (numel(nonzeros(cat(1,p0DetectedIndexVectors.L{:},p0DetectedIndexVectors.R{:})))/30)/numel(testSignal);
    evaluation.doaRatio = (numel(nonzeros(evaluation.azDeg))/30)/length(testSignal);
        
    for iBand = 1:30
        evaluation.targetp0DetectedIndexVectors.L{iBand} = p0DetectedIndexVectors.L{iBand}(targetSampleIndices.L{iBand});
        evaluation.interfp0DetectedIndexVectors.L{iBand} = p0DetectedIndexVectors.L{iBand}(interfSampleIndices.L{iBand});
        evaluation.targetp0DetectedIndexVectors.R{iBand} = p0DetectedIndexVectors.R{iBand}(targetSampleIndices.R{iBand});
        evaluation.interfp0DetectedIndexVectors.R{iBand} = p0DetectedIndexVectors.R{iBand}(interfSampleIndices.R{iBand});
    end
    for iBand = [3,4,5,10,14,20,25,30]
        figure;
        title('Gammatone band no.',num2str(iBand))

        subplot(1,2,1)
        title('Processed periodic samples - Left channel')
        xlabel('signal sample indices')
        ylabel('detected frequency (Hz)')
        hold on;
        scatter(targetSampleIndices.L{iBand}, evaluation.p0SearchRangeFreqVector(evaluation.targetp0DetectedIndexVectors.L{iBand}));
        scatter(interfSampleIndices.L{iBand}, evaluation.p0SearchRangeFreqVector(evaluation.interfp0DetectedIndexVectors.L{iBand}));
        legend('target','interferer')
        ylim(AlgorithmParameters.p0SearchRangeHz)
    
        subplot(1,2,2)
        title('Processed periodic samples - Right channel')
        xlabel('signal sample indices')
        ylabel('detected frequency (Hz)')
        hold on;
        scatter(targetSampleIndices.R{iBand}, evaluation.p0SearchRangeFreqVector(evaluation.targetp0DetectedIndexVectors.R{iBand}));
        scatter(interfSampleIndices.R{iBand}, evaluation.p0SearchRangeFreqVector(evaluation.interfp0DetectedIndexVectors.R{iBand}));
        legend('target','interferer')
        ylim(AlgorithmParameters.p0SearchRangeHz)
    end
end