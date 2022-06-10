% This is the MATLAB code for the master thesis "Speech enhacment for
% hearing aids using a combination of binaural and periodicity features" by
% Jeffrey Thomsen, written at the CvO University of Oldenburg, Germany, in
% 2022.
% This is the main script, the running of which will allow you to define
% custom simulation, algorithm and evaluation parameters
% The simulation is executed and the data and
% evaluation results are displayed and stored accordingly.

%% Preparation

clear

startTime = tic;

Publication = 'MA'; % 'DAGA', 'MA'

% Evaluation flags
Plotting = false;
BSIM = false;
Play = false;

% Initialize algorithm parameter and filter state structs
AlgorithmParameters = AlgorithmParametersConstructor();

% Alterations of simulation parameters
% main parameters of interest
AlgorithmParameters.snrCondition = false;
AlgorithmParameters.coherenceMask = false;

% to evaluate how much performance cancellation adds
AlgorithmParameters.Cancellation = true;
AlgorithmParameters.Enhancement = false;
AlgorithmParameters.RandomP0 = false;

% ideas for dealing with fc>1.4kHz
AlgorithmParameters.ChenP0Detection = false;
AlgorithmParameters.azimuthPooling = false;

% additional parameters
AlgorithmParameters.DOAProcessing = true;
AlgorithmParameters.bruemannSnrCondition = false;
AlgorithmParameters.bruemannSimpleDOA = false;

% Generate/load ITD/IPD-to-azimuth mapping function lookup table
tic
switch Publication
    case 'MA'
load('2022-05-26_itd_lookuptable_annotated.mat');
    case 'DAGA'
load('2022-03-04_itd_lookuptable_annotated.mat');
end
lookuptable = lookuptable.lookuptable;

% hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
% lookuptable = interauralToAzimuthLookuptable(hrtf, AlgorithmParameters,...
%     "0251M.flac","0652M.flac","1462F.flac","2035F.flac","2277F.flac",...
%     "3575F.flac","5694M.flac","7176M.flac","7729M.flac","7976F.flac");

AlgorithmParameters.lookuptable = lookuptable;
toc

% set sampling frequency to at least 2x upper gammatone centre frequency
samplingRateHz = ceil(2*AlgorithmParameters.Gammatone.fHigh);
AlgorithmParameters.Gammatone.samplingRateHz = samplingRateHz;

% generate gammatone filterbank
[AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);

centerFreqsHz = AlgorithmStates.L.GammatoneStates.analyzer.center_frequencies_hz;
nBands = AlgorithmParameters.Gammatone.nBands;
%%
switch Publication
    case 'MA'

%% Metadata for saving simulation data later

MetaData = struct;
[~,MetaData.gitCommitHash] = system('git rev-parse --short HEAD');
gitURL = 'https://github.com/jeffrey-thomsen/master-thesis-code/commit/';
MetaData.gitCommitURL = strcat(gitURL, MetaData.gitCommitHash);
MetaData.date = datetime('now','Format','yyyy-MM-dd''_''HH-mm');
dateString = string(MetaData.date);


%% Generate test signal battery

tic
% define parameters
TestSignalParameters.testSignalType = 'Battery';
TestSignalParameters.speakerIds = ...
    ["0251M", "0652M", "1462F", "2035F", "2277F", ...
     "3575F", "5694M", "7176M", "7729M", "7976F"];
TestSignalParameters.twoSpeakerVariety = true;
TestSignalParameters.targetAngles = [90, 0, -30];%[-90, -15, 0, 30, 60];
TestSignalParameters.nSpeakers = 2;

hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');

[inputMixedSignal, inputTargetSignal, inputInterfSignal, ...
    hrtfSamplingRateHz, inputMixedSignal_H, ...
    anglePermutations, speakerCombinations] = ...
    testSignalGenerator(TestSignalParameters, hrtf);

%% subset of only every second speakerCombination
inputMixedSignal = inputMixedSignal(1:2:end,:);
inputTargetSignal = inputTargetSignal(1:2:end,:);
inputInterfSignal = inputInterfSignal(1:2:end,:);
inputMixedSignal_H = inputMixedSignal_H(1:2:end,:);
speakerCombinations = speakerCombinations(1:2:end,:);

%%

% resample signals to algorithm sampling rate
for iSignal = 1:numel(inputMixedSignal)
    inputMixedSignal{iSignal}   = resample(inputMixedSignal{iSignal}, ...
        samplingRateHz, hrtfSamplingRateHz);
    inputMixedSignal_H{iSignal} = resample(inputMixedSignal_H{iSignal}, ...
        samplingRateHz, hrtfSamplingRateHz);
    inputTargetSignal{iSignal}  = resample(inputTargetSignal{iSignal}, ...
        samplingRateHz, hrtfSamplingRateHz);
    inputInterfSignal{iSignal}  = resample(inputInterfSignal{iSignal}, ...
        samplingRateHz, hrtfSamplingRateHz);
end

% time vector for plotting
dt = 1/AlgorithmParameters.Gammatone.samplingRateHz;
timeVec = dt:dt:dt*length(inputMixedSignal{1,1});
toc


%% Speech enhancement processing

nSignals = numel(inputMixedSignal);
nSpeakerCombos = size(speakerCombinations, 1);
nAnglePerms = size(anglePermutations, 1);

outputMixedSignal = cell(nSpeakerCombos, nAnglePerms);
SimulationData = cell(nSpeakerCombos, nAnglePerms);

estProcTime = 0;
elapsedTime = 0;
iSignal = 0;
for iSp = 1:nSpeakerCombos
    for jAn = 1:nAnglePerms

        iSignal = iSignal + 1;
        fprintf('Processing signal %1i of %1i\n\n', iSignal, nSignals);
        tic

        AlgorithmParameters.targetRangeDeg = ...
            anglePermutations(jAn, 1) + [-5 5];

        [outputMixedSignal{iSp, jAn}, ~, ...
            SimulationData{iSp, jAn}.Data] = ...
            speechEnhancement(inputMixedSignal{iSp, jAn}, ...
            AlgorithmParameters, AlgorithmStates);

        SimulationData{iSp, jAn}.signalIndices = ...
            1:length(inputMixedSignal{iSp, jAn});

        [estProcTime, elapsedTime] = estimateRemainingTime(iSignal, estProcTime, nSignals, elapsedTime);
    end
end


%% Evaluation using SimulationData (separated for memory reasons)

maskAppliedTarget = cell(nSpeakerCombos, nAnglePerms);
maskAppliedInterf = cell(nSpeakerCombos, nAnglePerms);
outputTargetSignal_S = cell(nSpeakerCombos, nAnglePerms);
outputInterfSignal_S = cell(nSpeakerCombos, nAnglePerms);

estProcTime = 0;
elapsedTime = 0;
iSignal = 0;
for iSp = 1:nSpeakerCombos
    for jAn = 1:nAnglePerms

        iSignal = iSignal + 1;
        fprintf('SimulationData Evaluating signal %1i of %1i\n\n', iSignal, nSignals);
        tic

        [maskAppliedTarget{iSp,jAn}, maskAppliedInterf{iSp,jAn}] = ...
            extractAppliedMask(SimulationData{iSp,jAn}, nBands);

        % compute separate enhanced signals for SNR improvement calculation
        [outputTargetSignal_S{iSp,jAn}, outputInterfSignal_S{iSp,jAn}] = ...
            shadowMethod(inputTargetSignal{iSp,jAn}, inputInterfSignal{iSp,jAn}, ...
            anglePermutations(jAn,1), SimulationData{iSp,jAn}.Data, ...
            AlgorithmParameters, AlgorithmStates);

        [estProcTime, elapsedTime] = estimateRemainingTime(iSignal, estProcTime, nSignals, elapsedTime);
    end
end


%% Save SimulationData subset

% tic
% 
% p0SearchRangeSamplesVector = SimulationData{1}.Data.p0SearchRangeSamplesVector;
% nSpeakerCombos = size(speakerCombinations, 1);
% nAnglePerms = size(anglePermutations, 1);
% SubsetSimulationData = cell(nSpeakerCombos, nAnglePerms);
% for iSp=1:nSpeakerCombos
%     for jAn=1:nAnglePerms
%         SubsetSimulationData{iSp, jAn}.p0DetectedIndexVectors = ...
%             SimulationData{iSp, jAn}.Data.p0DetectedIndexVectors;
%         SubsetSimulationData{iSp, jAn}.targetSampleIndices = ...
%             SimulationData{iSp, jAn}.Data.targetSampleIndices;
%         SubsetSimulationData{iSp, jAn}.interfSampleIndices = ...
%             SimulationData{iSp, jAn}.Data.interfSampleIndices;
%     end
% end
% filename = string(dateString+'_glimpse_data.mat');
% save(filename, 'SubsetSimulationData', 'MetaData', '-v7.3')
% 
% toc

tic

filename = string(dateString+'_simulation_data.mat');
save(filename, ...
'AlgorithmParameters', 'TestSignalParameters', ...
'anglePermutations', 'speakerCombinations', ...
'SimulationData', ...
'MetaData', '-v7.3')

toc


%%

clear SimulationData


%% Evaluation

% inputSnr = cell(nSpeakerCombos, nAnglePerms);
ibmTarget = cell(nSpeakerCombos, nAnglePerms);
ibmInterf = cell(nSpeakerCombos, nAnglePerms);

precision = cell(nSpeakerCombos, nAnglePerms);
recall = cell(nSpeakerCombos, nAnglePerms);
outputTargetSignal_H = cell(nSpeakerCombos, nAnglePerms);
outputInterfSignal_H = cell(nSpeakerCombos, nAnglePerms);

deltaSNR_H = zeros(nSpeakerCombos, nAnglePerms);
deltaSNR_S = zeros(nSpeakerCombos, nAnglePerms);
% outputSnr = cell(nSpeakerCombos, nAnglePerms);

estProcTime = 0;
elapsedTime = 0;
iSignal = 0;
for iSp = 1:nSpeakerCombos
    for jAn = 1:nAnglePerms

        iSignal = iSignal + 1;
        fprintf('Evaluating signal %1i of %1i\n\n', iSignal, nSignals);
        tic

        % compute IBM and extract glimpse masks from simulation data
        [inputTargetSubbandSignals, ~] = subbandDecompositionBinaural(...
            inputTargetSignal{iSp,jAn}, AlgorithmStates);
        [inputInterfSubbandSignals, ~] = subbandDecompositionBinaural(...
            inputInterfSignal{iSp,jAn}, AlgorithmStates);
        
        inputSnr.L = computeTimeFreqSnr(inputTargetSubbandSignals.L, ...
            inputInterfSubbandSignals.L, samplingRateHz, nBands);
        inputSnr.R = computeTimeFreqSnr(inputTargetSubbandSignals.R, ...
            inputInterfSubbandSignals.R, samplingRateHz, nBands);

        ibmTarget{iSp,jAn}.L = inputSnr.L >  0;
        ibmTarget{iSp,jAn}.R = inputSnr.R >  0;
        ibmInterf{iSp,jAn}.L = inputSnr.L <= 0;
        ibmInterf{iSp,jAn}.R = inputSnr.R <= 0;
        
        % compare chosen glimpses to IBM
        [precision{iSp,jAn}, recall{iSp,jAn}] = compareToIbm(...
            maskAppliedTarget{iSp,jAn}, maskAppliedInterf{iSp,jAn}, ...
            ibmTarget{iSp,jAn}, ibmInterf{iSp,jAn});

        if Plotting
            plotIbmGlimpses(maskAppliedTarget{iSp,jAn}, maskAppliedInterf{iSp,jAn}, ...
                ibmTarget{iSp,jAn}, ibmInterf{iSp,jAn}, inputTargetSignal{iSp,jAn}, ...
                inputInterfSignal{iSp,jAn}, SimulationData{iSp,jAn}, timeVec, ...
                AlgorithmParameters.Gammatone.samplingRateHz);
        end

        % compute separate enhanced signals for SNR improvement calculation
        [outputTargetSignal_H{iSp,jAn}, outputInterfSignal_H{iSp,jAn}] = ...
            hagermanMethod('pre-calc', anglePermutations(jAn,1), ...
            AlgorithmParameters, AlgorithmStates, ...
            outputMixedSignal{iSp,jAn}, inputMixedSignal_H{iSp,jAn});

        % compute SNR improvement
        deltaSNR_H(iSp,jAn) = computeSnrImprovement(...
            inputTargetSignal{iSp,jAn}, inputInterfSignal{iSp,jAn}, ...
            outputTargetSignal_H{iSp,jAn}, outputInterfSignal_H{iSp,jAn});
        deltaSNR_S(iSp,jAn) = computeSnrImprovement(...
            inputTargetSignal{iSp,jAn}, inputInterfSignal{iSp,jAn}, ...
            outputTargetSignal_S{iSp,jAn}, outputInterfSignal_S{iSp,jAn});

%         % compute TF gray mask SNR improvement
%         [outputTargetSubbandSignals_H, ~] = subbandDecompositionBinaural(...
%             outputTargetSignal_H{iSp,jAn}, AlgorithmStates);
%         [outputInterfSubbandSignals_H, ~] = subbandDecompositionBinaural(...
%             outputInterfSignal_H{iSp,jAn}, AlgorithmStates);
% 
%         outputSnr{iSp,jAn}.L = computeTimeFreqSnr(outputTargetSubbandSignals_H.L, ...
%             outputInterfSubbandSignals_H.L, samplingRateHz, nBands);
%         outputSnr{iSp,jAn}.R = computeTimeFreqSnr(outputTargetSubbandSignals_H.R, ...
%             outputInterfSubbandSignals_H.R, samplingRateHz, nBands);
        
        % estimate SII improvement with BSIM (optional)
        if BSIM
            [deltaSrt_H{iSp,jAn}, SrtIn_H{iSp,jAn}, SrtOut_H{iSp,jAn}] = ...
                computeSiiImprovement(...
                inputTargetSignal{iSp,jAn}, inputInterfSignal{iSp,jAn}, ...
                outputTargetSignal_H{iSp,jAn}, outputInterfSignal_H{iSp,jAn},...
                samplingRateHz, anglePermutations(jAn,2));
            
            [deltaSrt_S{iSp,jAn}, SrtIn_S{iSp,jAn}, SrtOut_S{iSp,jAn}] = ...
                computeSiiImprovement(...
                inputTargetSignal{iSp,jAn}, inputInterfSignal{iSp,jAn}, ...
                outputTargetSignal_S{iSp,jAn}, outputInterfSignal_S{iSp,jAn},...
                samplingRateHz, anglePermutations(jAn,2));
        end

        [estProcTime, elapsedTime] = estimateRemainingTime(iSignal, estProcTime, nSignals, elapsedTime);
    end
end


%% Evaluation evaluation

evaluationPlots(speakerCombinations, anglePermutations, precision, recall, deltaSNR_H);


%% Play signals

if Play
    box1 = msgbox('Play original signal (Ensure volume is adequately set)');
    waitfor(box1);
    soundsc(inputMixedSignal{iSp,jAn}, samplingRateHz);
    box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
        ' section of script (adjusting filter variables if necessary)']);
    waitfor(box2);
    soundsc(outputMixedSignal{iSp,jAn}, samplingRateHz);
end


%% Save simulation data

% audiowrite(['testInput','.wav'], mixedSignal, ...
%     AlgorithmParameters.Gammatone.samplingRateHz);
% audiowrite(['testJeffrey','.wav'], enhancedMixedSignal, ...
%     AlgorithmParameters.Gammatone.samplingRateHz);

tic
filename = string(dateString+'_signal_data.mat');
save(filename,...
'AlgorithmParameters', 'TestSignalParameters', 'AlgorithmStates',...
'anglePermutations', 'speakerCombinations',...
'centerFreqsHz', 'samplingRateHz', 'timeVec', ...
'inputMixedSignal', 'inputTargetSignal', 'inputInterfSignal', 'inputMixedSignal_H',...
'outputMixedSignal','outputTargetSignal_H', 'outputInterfSignal_H',...
'outputTargetSignal_S', 'outputInterfSignal_S',...
'lookuptable', 'MetaData',...
'-v7.3')
toc

tic
filename = string(dateString+'_evaluation_data.mat');
save(filename,...
'AlgorithmParameters', 'TestSignalParameters', ...
'anglePermutations', 'speakerCombinations',...
'ibmTarget', 'ibmInterf', ...
'maskAppliedTarget', 'maskAppliedInterf', ...
'precision', 'recall', 'deltaSNR_S', 'deltaSNR_H', ...
'MetaData', '-v7.3')
% 'inputSnr', 'outputSnr', 
toc

%%
    case 'DAGA'
%%
%% DAGA
%%
%% Generate test signal (DAGA)

% DAGA: target 0 degrees, interferer 60 degrees
TestSignalParameters.testSignalType = 'DAGA';
[inputMixedSignal, inputTargetSignal, inputInterfSignal, hrtfSamplingRateHz] = ...
    testSignalGenerator(TestSignalParameters);

% adjust sampling rate to gammatone filterbank
mixedSignal = resample(inputMixedSignal, ...
    AlgorithmParameters.Gammatone.samplingRateHz, hrtfSamplingRateHz);
inputTargetSignal = resample(inputTargetSignal, ...
    AlgorithmParameters.Gammatone.samplingRateHz, hrtfSamplingRateHz);
inputInterfSignal = resample(inputInterfSignal, ...
    AlgorithmParameters.Gammatone.samplingRateHz, hrtfSamplingRateHz);

% time vector for plotting
dt = 1/AlgorithmParameters.Gammatone.samplingRateHz;
timeVec = dt:dt:dt*length(mixedSignal);

%% Parameters
AlgorithmParameters.p0SearchRangeHz = [100 350];
AlgorithmParameters.snrThresholdInDb = 10;
AlgorithmParameters.targetRangeDeg = [-5 5];

%% Enhance mixed speech
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

%% Evaluate enhanced mixed speech
evalMixed = evaluateAlgorithm(dataMixed.p0DetectedIndexVectors, ...
    AlgorithmParameters, dataMixed.p0SearchRangeSamplesVector, ...
    mixedSignal, timeVec, dataMixed.ivsMaskCells, ...
    dataMixed.azimuthDegCells, dataMixed.targetSampleIndices, ...
    dataMixed.interfSampleIndices, centerFreqsHz, false);

%% Enhance target speech
tic
[enhancedSignalTarget, AlgorithmStates, dataTarget] = ...
    speechEnhancement(inputTargetSignal, AlgorithmParameters, AlgorithmStates);
toc
enhancedSignalTarget = enhancedSignalTarget./max(max(abs(enhancedSignalTarget)));

%% Evaluate enhanced target speech
evalTarget = evaluateAlgorithm(dataTarget.p0DetectedIndexVectors, ...
    AlgorithmParameters, dataTarget.p0SearchRangeSamplesVector, ...
    inputTargetSignal, timeVec, dataTarget.ivsMaskCells, ...
    dataTarget.azimuthDegCells, dataTarget.targetSampleIndices, ...
    dataTarget.interfSampleIndices, centerFreqsHz, false);

p0SearchRangeHz = evalMixed.p0SearchRangeFreqVector;

%% Enhance interferer speech
tic
[enhancedSignalInterf, AlgorithmStatesInterf, dataInterf] = ...
    speechEnhancement(inputInterfSignal, AlgorithmParameters, AlgorithmStates);
toc
enhancedSignalInterf = enhancedSignalInterf./max(max(abs(enhancedSignalInterf)));

%% Evaluate enhanced interferer speech
evalInterf = evaluateAlgorithm(dataInterf.p0DetectedIndexVectors, ...
    AlgorithmParameters, dataInterf.p0SearchRangeSamplesVector, ...
    inputInterfSignal, timeVec, dataInterf.ivsMaskCells, ...
    dataInterf.azimuthDegCells, dataInterf.targetSampleIndices, ...
    dataInterf.interfSampleIndices, centerFreqsHz, false);

%% DAGA results

[precision, recall] = glimpseDistanceMetric(AlgorithmParameters, ...
    dataMixed, dataTarget, dataInterf, centerFreqsHz, mixedSignal);

plotDAGAresults(evalMixed, evalTarget, evalInterf, ...
    mixedSignal, inputTargetSignal, inputInterfSignal, ...
    dataMixed, dataTarget, dataInterf, AlgorithmParameters, timeVec)

%% Play mixed signals
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
soundsc(mixedSignal, AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
soundsc(enhancedSignalMixed, AlgorithmParameters.Gammatone.samplingRateHz);

%% Play target signals
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
soundsc(inputTargetSignal, AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
soundsc(enhancedSignalTarget, AlgorithmParameters.Gammatone.samplingRateHz);

%% Play interferer signals
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
soundsc(inputInterfSignal, AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
soundsc(enhancedSignalInterf, AlgorithmParameters.Gammatone.samplingRateHz);

%% Save simulation data

% audiowrite(['testInput','.wav'], mixedSignal, ...
%     AlgorithmParameters.Gammatone.samplingRateHz);
% audiowrite(['testJeffrey','.wav'], enhancedMixedSignal, ...
%     AlgorithmParameters.Gammatone.samplingRateHz);

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
'inputInterfSignal', 'mixedSignal', 'inputTargetSignal',...
'lookuptable',...
'MetaData')

end


%%

endTime = toc(startTime);
fprintf('Total elapsed time (hh:mm:ss): %02.0f:%02.0f:%02.0f \n', ...
            floor(endTime/3600), ...
            floor(mod(endTime,3600)/60), ...
            floor(mod(endTime,60)));


%% Auxiliary evaluation functions

function [estProcTime, elapsedTime] = ...
    estimateRemainingTime(iSignal, estProcTime, nSignals, elapsedTime)
    
    if iSignal > 1
        estProcTime = 0.9*estProcTime + 0.1*toc;
    else
        estProcTime = toc;
    end
    estRemTime1 = (nSignals-iSignal)*estProcTime;

    elapsedTime = elapsedTime + toc;
    estRemTime2 = (nSignals-iSignal)*(elapsedTime/iSignal);
    
    fprintf('Estimated remaining time (hh:mm:ss): %02.0f:%02.0f:%02.0f or %02.0f:%02.0f:%02.0f\n', ...
        floor(estRemTime1/3600), ...
        floor(mod(estRemTime1,3600)/60), ...
        floor(mod(estRemTime1,60)), ...
        floor(estRemTime2/3600), ...
        floor(mod(estRemTime2,3600)/60), ...
        floor(mod(estRemTime2,60)));
end

function evaluation = evaluateAlgorithm(p0DetectedIndexVectors, ...
    AlgorithmParameters, p0SearchRangeSamplesVector, testSignal, timeVec, ...
    ivsMask, azimuthDegCells, targetSampleIndices, interfSampleIndices, fc,...
    Plotting)

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
     nnz(evaluation.coherentPeriodicComponentLogicalVector.R)) /...
     (2*numel(evaluation.coherentPeriodicComponentLogicalVector.L));

    %% p0 detection histogram
    [evaluation.GC,evaluation.GR] = groupcounts(cat(1,p0DetectedIndexVectors.L{:},...
        p0DetectedIndexVectors.R{:}));
    evaluation.GR = nonzeros(evaluation.GR);%+p0SearchRangeSamplesVector(1)-1;
    evaluation.GC = evaluation.GC(end-length(evaluation.GR)+1:end);
    evaluation.p0SearchRangeSecondsVector = ...
        p0SearchRangeSamplesVector./AlgorithmParameters.Gammatone.samplingRateHz;
    evaluation.p0SearchRangeFreqVector = ...
        AlgorithmParameters.Gammatone.samplingRateHz./p0SearchRangeSamplesVector;
    
    if Plotting
    figure;
    hBar = bar(evaluation.p0SearchRangeSecondsVector(evaluation.GR),...
        evaluation.GC);
%     set(gca,'XScale','log');
    title('p0 detection histogram',['total no. of subband samples: ',...
        num2str(numel(testSignal)*nBands)])
    xlabel('Period (s)')
    ylabel('No. of occurences in subband samples')
    set(gca,'FontSize',14);
    end
    %% signal spectrogram
    if Plotting
    
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
    end

    %% plot azimuth estimation histogram
    if Plotting
    figure;
    histogram(nonzeros(evaluation.azDeg(abs(evaluation.azDeg)<90)),...
        'Normalization','probability','NumBins',90)
    title(['DOA estimation histogram'])%,['total no. of subband samples: ',num2str(length(testSignal)*nBands)])
%     ylabel({'No. of occurences';'in subband samples'})
    ylabel('Probability of occurrence')
    xlabel('Azimuth (°)')
    xticks([-90 -60 -30 0 30 60 90])
%     ylim([0 0.25])
    set(gca,'FontSize',14);
    end
    %% plot position and p0 of detected target and interferer samples
%     iStart = 1;
%     iEnd = 30;
%     if Plotting
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
%     end
end

