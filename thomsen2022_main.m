%% Start/Simulation interface
clear all;
msgbox('Hi, you are running the Thomsen2022 speech enhancement algorithm!')

TestSignalParameters = struct;
BlockFeedingParameters = struct;
AlgorithmParameters = struct;
%% Test signal generator

testSignal = testSignalGenerator;
testSignal = testInputSignal(testSignal);
%% Block-feeding routine

processedSignal = blockFeedingRoutine(testSignal);
