%% Main function to generate tests
function tests = basicTest
tests = functiontests(localfunctions);
end

%% Test functions

% function testMainScript(testCase)
%     thomsen2022_main;
% end

function testDefaultTestSignalGeneratorSignalOutputLength(testCase)
    testSignal = testSignalGenerator;
    verifyEqual(testCase, length(testSignal), 1001, "RelTol", 0.01)
    verifyLessThanOrEqual(testCase, length(testSignal), 1000)
    expectedSize = [1000 2];
    verifySize(testCase,testSignal,expectedSize)
end

function testPredefinedTestSignalGeneratorSignalOutputLength(testCase)
    nSamples = 3000;
    TestSignalParameters.nSamples = nSamples;
    testSignal = testSignalGenerator(TestSignalParameters);
    expectedSize = [nSamples 2];
    verifySize(testCase,testSignal,expectedSize)
end

function testTestInputSignalConvertsArrayDimensions(testCase)
    input = rand(2,10000)-0.5;
    output = testInputSignal(input);
    expectedSize = [10000 2];
    verifySize(testCase,output,expectedSize)
end

function testTestInputSignalOutputEqualsInput(testCase)
    input = 2*rand(10000,2)-1;
    output = testInputSignal(input);
    verifyEqual(testCase,output,input)
end

function testSpeechEnhancementOutputEqualsInput(testCase)
    input = rand(10000,2)-0.5;
    output = speechEnhancement(input);
    verifyEqual(testCase,output,input)
end

function testDefaultBlockFeedingRoutineOutputEqualsInput(testCase)
    input = rand(10000,2)-0.5;
    output = blockFeedingRoutine(input);
    verifyEqual(testCase,output,input)
end

function testBlockFeedingRoutineOutputEqualsInput(testCase)
    input = rand(10000,2)-0.5;
    BlockFeedingParameters.blockLength = 100;
    output = blockFeedingRoutine(input, BlockFeedingParameters);
    verifyEqual(testCase,output,input)
end

function testBlockFeedingRoutineOutputEqualsInput1SampleAtATime(testCase)
    input = rand(10000,2)-0.5;
    BlockFeedingParameters.blockLength = 1;
    output = blockFeedingRoutine(input, BlockFeedingParameters);
    verifyEqual(testCase,output,input)
end