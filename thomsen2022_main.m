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
msgbox('Hi, you are running the Thomsen2022 speech enhancement algorithm!')

%% Preparation
TestSignalParameters = struct;
TargetAngleParameters = struct;
BlockFeedingParameters = struct;
AlgorithmParameters = AlgorithmParametersConstructor();

hrtf = SOFAload('HRIR_KEMAR_DV0001_3.sofa',[5 2],'R');
AlgorithmParameters.Gammatone.samplingRateHz = hrtf.Data.SamplingRate;
AlgorithmParameters.lookuptable = ...
    ipdToAzimuthLookuptable(hrtf, AlgorithmParameters);

[AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);


%% Test signal generator
testSignal = testSignalGenerator;
testSignal = testInputSignal(testSignal);

%% Block-feeding routine
processedSignal = blockFeedingRoutine(testSignal,BlockFeedingParameters,...
    AlgorithmParameters, AlgorithmStates);

%% Evaluation

%% Display and save data
