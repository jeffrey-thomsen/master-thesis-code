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

% load HRTF
hrtf = SOFAload('HRIR_KEMAR_DV0001_3.sofa',[5 2],'R');
AlgorithmParameters.Gammatone.samplingRateHz = hrtf.Data.SamplingRate;

% generate IPD-to-azimuth mapping function
AlgorithmParameters.lookuptable = ...
    ipdToAzimuthLookuptable(hrtf, AlgorithmParameters);

% resampling sampling frequency to at least 
% 2x upper gammatone centre frequency
resampleFactor = floor(AlgorithmParameters.Gammatone.samplingRateHz/...
    (2*AlgorithmParameters.Gammatone.fHigh));
AlgorithmParameters.Gammatone.samplingRateHz =...
    AlgorithmParameters.Gammatone.samplingRateHz / resampleFactor; 

% generate gammatone filterbank
[AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);

%% Generate test signal

% load clean speech signals
[target,fs]=audioread('sp01.wav');
target = resample(target,hrtf.Data.SamplingRate,fs);

[interferer,fs]=audioread('sp30.wav');
interferer = resample(interferer,hrtf.Data.SamplingRate,fs);

% equalize lengths
target = target(1:length(interferer)); 

% equalize levels - SNR 0dB
interferer = interferer*(sum(abs(target).^2)./sum(abs(interferer).^2));

% convolve with HRTF at different incidence angles
target = SOFAspat(target, hrtf, 0, 0);
interferer = SOFAspat(interferer, hrtf, 60, 0);

% add, resample and normalize test signal
testSignal = target + interferer;
testSignal = resample(testSignal,1,resampleFactor); 
testSignal = testSignal./max(max(testSignal));
% testSignal = testSignalGenerator;
testSignal = testInputSignal(testSignal);

%% Block-feeding routine - Speech enhancement algorithm
% tic
% processedSignal = blockFeedingRoutine(testSignal,BlockFeedingParameters,...
%     AlgorithmParameters, AlgorithmStates);
% toc
tic
[enhancedSignal, AlgorithmStates, p0DetectedIndexVectors, ...
    p0SearchRangeSamplesVector, ipdRad, ivsMask, ...
    ipdDisambiguatedLogicalCells, azimuthDegCells,...
    targetSampleIndices, interfererSampleIndices] = ...
    speechEnhancement(testSignal, AlgorithmParameters, AlgorithmStates);
toc
%% Evaluation

[GC,GR] = groupcounts(cat(1,p0DetectedIndexVectors.L{:},p0DetectedIndexVectors.R{:}));
GR = nonzeros(GR);%+p0SearchRangeSamplesVector(1)-1;
GC = GC(end-length(GR)+1:end);
p0SearchRangeFreqVector = AlgorithmParameters.Gammatone.samplingRateHz./p0SearchRangeSamplesVector;
% linspace(AlgorithmParameters.p0SearchRangeHz(1),AlgorithmParameters.p0SearchRangeHz(2),length(p0SearchRangeSamplesVector));
figure;hBar=bar(p0SearchRangeFreqVector(GR),GC);set(gca,'XScale','log');
meanTestSignal = (testSignal(:,1) + testSignal(:,2))./2;
figure;spectrogram(meanTestSignal,hamming(1000),[],[],AlgorithmParameters.Gammatone.samplingRateHz);xlim([0, 0.5]);

azDeg = zeros(30,length(ivsMask{1}));
for iBand = 1:30
    azDeg(iBand,ivsMask{iBand}) = azimuthDegCells{iBand};
    coherentPeriodicComponentLogicalVectorL(iBand,:) = ...
        ivsMask{iBand} & p0DetectedIndexVectors.L{iBand};
    coherentPeriodicComponentLogicalVectorR(iBand,:) = ...
        ivsMask{iBand} & p0DetectedIndexVectors.R{iBand};
end
figure;histogram(nonzeros(azDeg(abs(azDeg)<90)))

periodicityRatio = (numel(nonzeros(cat(1,p0DetectedIndexVectors.L{:},p0DetectedIndexVectors.R{:})))/30)/numel(testSignal);
doaRatio = (numel(nonzeros(azDeg))/30)/length(testSignal);


%% Compare to target speech
targetSignal = resample(target,1,resampleFactor);
tic
[enhancedTargetSignal, AlgorithmStates, p0DetectedIndexVectors, ...
    p0SearchRangeSamplesVector, ipdRad, ivsMask, ...
    ipdDisambiguatedLogicalCells, azimuthDegCells...
    , targetSampleIndicesTar, interfererSampleIndicesTar] = ...
    speechEnhancement(targetSignal, AlgorithmParameters, AlgorithmStates);
toc
% Evaluation

[GCTar,GRTar] = groupcounts(cat(1,p0DetectedIndexVectors.L{:},p0DetectedIndexVectors.R{:}));
GRTar = nonzeros(GRTar);%+p0SearchRangeSamplesVector(1)-1;
GCTar = GCTar(end-length(GRTar)+1:end);
p0SearchRangeFreqVector = AlgorithmParameters.Gammatone.samplingRateHz./p0SearchRangeSamplesVector;
% linspace(AlgorithmParameters.p0SearchRangeHz(1),AlgorithmParameters.p0SearchRangeHz(2),length(p0SearchRangeSamplesVector));
figure;hBar=bar(p0SearchRangeFreqVector(GRTar),GCTar);set(gca,'XScale','log');
meanTargetSignal = (targetSignal(:,1) + targetSignal(:,2))./2;
figure;spectrogram(meanTargetSignal,hamming(1000),[],[],AlgorithmParameters.Gammatone.samplingRateHz);xlim([0, 0.5]);

azDegTar = zeros(30,length(ivsMask{1}));
for iBand = 1:30
    azDegTar(iBand,ivsMask{iBand}) = azimuthDegCells{iBand};
    coherentPeriodicComponentLogicalVectorLTarget(iBand,:) = ...
        ivsMask{iBand} & p0DetectedIndexVectors.L{iBand};
    coherentPeriodicComponentLogicalVectorRTarget(iBand,:) = ...
        ivsMask{iBand} & p0DetectedIndexVectors.R{iBand};
end
figure;histogram(nonzeros(azDegTar(abs(azDegTar)<90)))

periodicityRatioTarget = (numel(nonzeros(cat(1,p0DetectedIndexVectors.L{:},p0DetectedIndexVectors.R{:})))/30)/numel(targetSignal);
doaRatioTarget = (numel(nonzeros(azDegTar))/30)/length(targetSignal);


%% Compare to interferer speech
interfSignal = resample(interferer,1,resampleFactor);
tic
[enhancedInterfSignal, AlgorithmStates, p0DetectedIndexVectors, ...
    p0SearchRangeSamplesVector, ipdRad, ivsMask, ...
    ipdDisambiguatedLogicalCells, azimuthDegCells...
    , targetSampleIndicesInt, interfererSampleIndicesInt] = ...
    speechEnhancement(interfSignal, AlgorithmParameters, AlgorithmStates);
toc
% Evaluation

[GCInt,GRInt] = groupcounts(cat(1,p0DetectedIndexVectors.L{:},p0DetectedIndexVectors.R{:}));
GRInt = nonzeros(GRInt);%+p0SearchRangeSamplesVector(1)-1;
GCInt = GCInt(end-length(GRInt)+1:end);
p0SearchRangeFreqVector = AlgorithmParameters.Gammatone.samplingRateHz./p0SearchRangeSamplesVector;
% linspace(AlgorithmParameters.p0SearchRangeHz(1),AlgorithmParameters.p0SearchRangeHz(2),length(p0SearchRangeSamplesVector));
figure;hBar=bar(p0SearchRangeFreqVector(GRInt),GCInt);set(gca,'XScale','log');
meanInterfSignal = (interfSignal(:,1) + interfSignal(:,2))./2;
figure;spectrogram(meanInterfSignal,hamming(1000),[],[],AlgorithmParameters.Gammatone.samplingRateHz);xlim([0, 0.5]);

azDegInt = zeros(30,length(ivsMask{1}));
for iBand = 1:30
    azDegInt(iBand,ivsMask{iBand}) = azimuthDegCells{iBand};
    coherentPeriodicComponentLogicalVectorLInterf(iBand,:) = ...
        ivsMask{iBand} & p0DetectedIndexVectors.L{iBand};
    coherentPeriodicComponentLogicalVectorRInterf(iBand,:) = ...
        ivsMask{iBand} & p0DetectedIndexVectors.R{iBand};
end
figure;histogram(nonzeros(azDegInt(abs(azDegInt)<90)))

periodicityRatioInterf = (numel(nonzeros(cat(1,p0DetectedIndexVectors.L{:},p0DetectedIndexVectors.R{:})))/30)/numel(interfSignal);
doaRatioInterf = (numel(nonzeros(azDegInt))/30)/length(interfSignal);

%% Compare coherent periodic samples detected

% overlap of target and interferer detected samples in mixed signal
overlap=(coherentPeriodicComponentLogicalVectorL&coherentPeriodicComponentLogicalVectorLTarget)...% detected samples detected in target and mixed signal
&...
(coherentPeriodicComponentLogicalVectorL&coherentPeriodicComponentLogicalVectorLInterf); % detected samples detected in interferer and mixed signal

% detected samples in mixed signal that are neither in target nor in
% interferer signal
nooverlap=(coherentPeriodicComponentLogicalVectorL&~coherentPeriodicComponentLogicalVectorLTarget)...
&...
(coherentPeriodicComponentLogicalVectorL&~coherentPeriodicComponentLogicalVectorLInterf);

% detected samples in mixed signal that are either in target or interferer
% signal
effectiveDetection = (coherentPeriodicComponentLogicalVectorL&coherentPeriodicComponentLogicalVectorLTarget)...% detected samples detected in target and mixed signal
|...
(coherentPeriodicComponentLogicalVectorL&coherentPeriodicComponentLogicalVectorLInterf); % detected samples detected in interferer and mixed signal


proportionTargetInMixed = numel(nonzeros(coherentPeriodicComponentLogicalVectorL&coherentPeriodicComponentLogicalVectorLTarget))/numel(nonzeros(coherentPeriodicComponentLogicalVectorL));
proportionInterfInMixed = numel(nonzeros(coherentPeriodicComponentLogicalVectorL&coherentPeriodicComponentLogicalVectorLInterf))/numel(nonzeros(coherentPeriodicComponentLogicalVectorL));
proportionTargetInterfOverlapInMixed = numel(nonzeros(overlap))/numel(nonzeros(coherentPeriodicComponentLogicalVectorL));



% proportionMixedInTarget = numel(nonzeros(coherentPeriodicComponentLogicalVectorL&coherentPeriodicComponentLogicalVectorLTarget))/numel(nonzeros(coherentPeriodicComponentLogicalVectorLTarget));

% mixed processed samples that were not processed in target or interferer
proportionMixedNotInTargetOrInterf = numel(nonzeros(nooverlap))/numel(nonzeros(coherentPeriodicComponentLogicalVectorL));
% complementary value:
proportionTargetInterfInMixed = numel(nonzeros(effectiveDetection))/numel(nonzeros(coherentPeriodicComponentLogicalVectorL));

%
tSI = numel(cat(1,targetSampleIndices.L{:},targetSampleIndices.R{:}));
tSIT = numel(cat(1,targetSampleIndicesTar.L{:},targetSampleIndicesTar.R{:}));
tSII = numel(cat(1,targetSampleIndicesInt.L{:},targetSampleIndicesInt.R{:}));

iSI = numel(cat(1,interfererSampleIndices.L{:},interfererSampleIndices.R{:}));
iSIT = numel(cat(1,interfererSampleIndicesTar.L{:},interfererSampleIndicesTar.R{:}));
iSII = numel(cat(1,interfererSampleIndicesInt.L{:},interfererSampleIndicesInt.R{:}));
%%
% Play signals
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
sound(testSignal,AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
sound(enhancedSignal,AlgorithmParameters.Gammatone.samplingRateHz);
%% Display and save data
% audiowrite(['testInput','.wav'],testSignal,AlgorithmParameters.Gammatone.samplingRateHz);
% audiowrite(['testJeffrey','.wav'],processedSignal,AlgorithmParameters.Gammatone.samplingRateHz);