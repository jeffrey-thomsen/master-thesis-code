%% Main function to generate tests
function tests = testSuite
tests = functiontests(localfunctions);
end

%% Test functions

% function testMainScript(testCase)
%     thomsen2022_main;
% end

function testDefaultTestSignalGeneratorSignalOutputLength(testCase)
% test if testSignalGenerator.m returns a signal of the expected size as
% default
    
    % Exercise
    testSignal = testSignalGenerator;
    % Verify
    verifyEqual(testCase, length(testSignal), 32001, "RelTol", 0.01)
    verifyLessThanOrEqual(testCase, length(testSignal), 32000)
    expectedSize = [32000 2];
    verifySize(testCase,testSignal,expectedSize)
end

function testPredefinedTestSignalGeneratorSignalOutputLength(testCase)
% test if testSignalGenerator.m returns a signal of specified length
    
    % Setup
    nSamples = 3000;
    TestSignalParameters.nSamples = nSamples;
    % Exercise
    testSignal = testSignalGenerator(TestSignalParameters);
    % Verify
    expectedSize = [nSamples 2];
    verifySize(testCase,testSignal,expectedSize)
end

function testTestInputSignalConvertsArrayDimensions(testCase)
% test if testInputSignal.m returns signal as Nx2 matrix
    
    % Setup
    input = rand(2,10000)-0.5;
    % Exercise
    output = testInputSignal(input);
    % Verify
    expectedSize = [10000 2];
    verifySize(testCase,output,expectedSize)
end

function testTestInputSignalOutputEqualsInput(testCase)
% test if testInputSignal.m leaves signal untouched if dimensions correct
    
    % Setup
    input = 2*rand(10000,2)-1;
    % Exercise
    output = testInputSignal(input);
    % Verify
    verifyEqual(testCase,output,input)
end

function testConstructGammatoneFilterbankExecutes(testCase)
% test if constructGammatoneFilterbank.m correctly outputs 2 structs
    
    % Setup
    Parameters = AlgorithmParametersConstructor();
    % Exercise
    [analyzer, synthesizer] = constructGammatoneFilterbank(Parameters.GammatoneParameters);
    % Verify
    verifyClass(testCase,analyzer,"struct")
    verifyClass(testCase,synthesizer,"struct")
end

function testGammatoneFilterbank(testCase)
% test if gammatone filterbank does near-perfect reconstruction of signal
    
    % Setup
    Parameters = AlgorithmParametersConstructor();
    % extract expected delay
    delaySeconds = Parameters.GammatoneParameters.desiredDelayInSeconds;
    samplingrateHz = Parameters.GammatoneParameters.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    [analyzer, synthesizer] = constructGammatoneFilterbank(Parameters.GammatoneParameters);
    % create test signal
    input = rand(10000,1)-0.5;

    % Exercise
    [decomposedSubbandSignals, ~] = ...
        subbandDecomposition(input, analyzer);
    [output, ~] = ...
        subbandResynthesis(decomposedSubbandSignals, synthesizer);
    
    % Verify
    % actual delay should correspond to desired delay
    actualDelay = finddelay(input,output);
    verifyEqual(testCase,actualDelay,delaySamples);
    
    % largest cross-correlation coefficient should be significantly larger
    % than second-largest
    sortedxCorr = sort(xcorr(input,output));
    xcorrRatio = sortedxCorr(end)/sortedxCorr(end-1);
    verifyGreaterThan(testCase,xcorrRatio,9.5)
end

function testGammatoneFilterbankInAlgorithm(testCase)
% test if filterbank does near-perfect reconstruction of signal within the
% algorithm framework
    
    % Setup
    AlgorithmParameters = AlgorithmParametersConstructor();
    % disable processing stage
    AlgorithmParameters.Cancellation = false;
    AlgorithmParameters.Enhancement = false;
    % extract expected delay
    delaySeconds = ...
        AlgorithmParameters.GammatoneParameters.desiredDelayInSeconds;
    samplingrateHz = ...
        AlgorithmParameters.GammatoneParameters.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    % create test signal
    input = rand(10000,2)-0.5;
    input = testInputSignal(input);
    
    % Exercise
    [output, ~] = ...
        speechEnhancement(input, AlgorithmParameters);
    
    % Verify
    for ch=1:2
        % actual delay should correspond to desired delay
        actualDelay = finddelay(input,output);
        verifyEqual(testCase,actualDelay(ch),delaySamples);

        % largest cross-correlation coefficient should be significantly 
        % larger than second-largest
        sortedxCorr = sort(xcorr(input(:,ch),output(:,ch)));
        xcorrRatio = sortedxCorr(end)/sortedxCorr(end-1);
        verifyGreaterThan(testCase,xcorrRatio,9.5)
    end
end

function testGammatoneFilterbankInBlockFeederDefault(testCase)
% test if filterbank does near-perfect reconstruction of signal within the
% block feeding framework with default block length
    
    % Setup
    BlockFeedingParameters = struct;
    AlgorithmParameters = AlgorithmParametersConstructor();
    % disable processing stage
    AlgorithmParameters.Cancellation = false;
    AlgorithmParameters.Enhancement = false;
    % extract expected delay
    delaySeconds = ...
        AlgorithmParameters.GammatoneParameters.desiredDelayInSeconds;
    samplingrateHz = ...
        AlgorithmParameters.GammatoneParameters.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    % create test signal
    input = rand(10000,2)-0.5;
    input = testInputSignal(input);

    % Exercise
    output = blockFeedingRoutine(input,BlockFeedingParameters,...
        AlgorithmParameters);
    
    % Verify
    for ch=1:2
        % actual delay should correspond to desired delay
        actualDelay = finddelay(input,output);
        verifyEqual(testCase,actualDelay(ch),delaySamples);

        % largest cross-correlation coefficient should be significantly 
        % larger than second-largest
        sortedxCorr = sort(xcorr(input(:,ch),output(:,ch)));
        xcorrRatio = sortedxCorr(end)/sortedxCorr(end-1);
        verifyGreaterThan(testCase,xcorrRatio,9.5)
    end
end

function testGammatoneFilterbankInBlockFeeder(testCase)
% test if filterbank does near-perfect reconstruction of signal within the
% block feeding framework with custom block length
    
    % Setup
    BlockFeedingParameters.blockLength = 2000;
    AlgorithmParameters = AlgorithmParametersConstructor();
    % disable processing stage
    AlgorithmParameters.Cancellation = false;
    AlgorithmParameters.Enhancement = false;
    % extract expected delay
    delaySeconds = ...
        AlgorithmParameters.GammatoneParameters.desiredDelayInSeconds;
    samplingrateHz = ...
        AlgorithmParameters.GammatoneParameters.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    % create test signal
    input = rand(10000,2)-0.5;
    input = testInputSignal(input);

    % Exercise
    output = blockFeedingRoutine(input,BlockFeedingParameters,...
        AlgorithmParameters);
    
    % Verify
    for ch=1:2
        % actual delay should correspond to desired delay
        actualDelay = finddelay(input,output);
        verifyEqual(testCase,actualDelay(ch),delaySamples);

        % largest cross-correlation coefficient should be significantly 
        % larger than second-largest
        sortedxCorr = sort(xcorr(input(:,ch),output(:,ch)));
        xcorrRatio = sortedxCorr(end)/sortedxCorr(end-1);
        verifyGreaterThan(testCase,xcorrRatio,9.5)
    end
end

function testGammatoneFilterbankInBlockFeeder1SampleAtATime(testCase)
% test if filterbank does near-perfect reconstruction of signal within the
% block feeding framework, running only one sample at a time
    
    % Setup
    BlockFeedingParameters.blockLength = 1;
    AlgorithmParameters = AlgorithmParametersConstructor();
    % disable processing stage
    AlgorithmParameters.Cancellation = false;
    AlgorithmParameters.Enhancement = false;
    % extract expected delay
    delaySeconds = ...
        AlgorithmParameters.GammatoneParameters.desiredDelayInSeconds;
    samplingrateHz = ...
        AlgorithmParameters.GammatoneParameters.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    % create test signal
    input = rand(10000,2)-0.5;
    input = testInputSignal(input);

    % Exercise
    output = blockFeedingRoutine(input,BlockFeedingParameters,...
        AlgorithmParameters);
    
    % Verify
    for ch=1:2
        % actual delay should correspond to desired delay
        actualDelay = finddelay(input,output);
        verifyEqual(testCase,actualDelay(ch),delaySamples);

        % largest cross-correlation coefficient should be significantly 
        % larger than second-largest
        sortedxCorr = sort(xcorr(input(:,ch),output(:,ch)));
        xcorrRatio = sortedxCorr(end)/sortedxCorr(end-1);
        verifyGreaterThan(testCase,xcorrRatio,9.5)
    end
end