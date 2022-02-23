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
tic
processedSignal = blockFeedingRoutine(testSignal,BlockFeedingParameters,...
    AlgorithmParameters, AlgorithmStates);
toc
%% Evaluation
% Play signals
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
sound(testSignal,AlgorithmParameters.Gammatone.samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
sound(processedSignal,AlgorithmParameters.Gammatone.samplingRateHz);
%% Display and save data
% audiowrite(['testInput','.wav'],testSignal,AlgorithmParameters.Gammatone.samplingRateHz);
% audiowrite(['testJeffrey','.wav'],processedSignal,AlgorithmParameters.Gammatone.samplingRateHz);