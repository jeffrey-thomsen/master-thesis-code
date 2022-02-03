%% Start/Simulation interface
clear
msgbox('Hi, you are running the Thomsen2022 speech enhancement algorithm!')

%% Preparation
TestSignalParameters = struct;
TargetAngleParameters = struct;
BlockFeedingParameters = struct;
AlgorithmParameters = AlgorithmParametersConstructor();

%% Test signal generator
testSignal = testSignalGenerator;
testSignal = testInputSignal(testSignal);

%% Block-feeding routine
processedSignal = blockFeedingRoutine(testSignal,BlockFeedingParameters,...
    AlgorithmParameters);
