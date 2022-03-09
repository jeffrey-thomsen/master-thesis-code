% These are the unit tests used for the development of the Thomsen2022
% speech enhancement algorithm. The test suite is run by executing
% "runtests" in MATLAB in the project folder.

%% Main function to generate tests
function tests = testSuite
tests = functiontests(localfunctions);
end

%% Test functions

% function testMainScript(testCase)
%     thomsen2022_main;
% end

% testSignalGenerator.m
function testDefaultTestSignalGeneratorSignalOutputLength(testCase)
% test if testSignalGenerator.m returns a signal of the expected size as
% default
    
    % Exercise
    testSignal = testSignalGenerator;
    % Verify
    verifyEqual(testCase, size(testSignal,1), 32001, "RelTol", 0.01)
    verifyLessThanOrEqual(testCase, size(testSignal,1), 32000)
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

% testInputSignal.m
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

% Gammatone Filterbank
function testConstructGammatoneFilterbankExecutes(testCase)
% test if constructGammatoneFilterbank.m correctly outputs 2 structs
% to contain the fiter parameters of the analysis and synthesis filterbank
    
    % Setup
    Parameters = AlgorithmParametersConstructor();
    % Exercise
    [analyzer, synthesizer] = constructGammatoneFilterbank(Parameters.Gammatone);
    % Verify
    verifyClass(testCase,analyzer,"struct")
    verifyClass(testCase,synthesizer,"struct")
end

function testGammatoneFilterbank(testCase)
% test if gammatone filterbank does near-perfect reconstruction of signal
    
    % Setup
    Parameters = AlgorithmParametersConstructor();
    % extract expected delay
    delaySeconds = Parameters.Gammatone.desiredDelayInSeconds;
    samplingrateHz = Parameters.Gammatone.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    [analyzer, synthesizer] = constructGammatoneFilterbank(Parameters.Gammatone);
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

function testGammatoneFilterbankBinaural(testCase)
% test if gammatone filterbank does near-perfect reconstruction of signal
% in the binaural framework
    
    % Setup
    Parameters = AlgorithmParametersConstructor();
    [States, Parameters.Gammatone.nBands] = ...
        AlgorithmStatesConstructor(Parameters);

    % extract expected delay
    delaySeconds = Parameters.Gammatone.desiredDelayInSeconds;
    samplingrateHz = Parameters.Gammatone.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    delaySamples = [delaySamples delaySamples];

    % create test signal
    input = testSignalGenerator;

    % Exercise
    [subbandSignalArray, ~] = ...
        subbandDecompositionBinaural(input, States);
    [output, ~] = ...
        subbandResynthesisBinaural(subbandSignalArray, States);
    
    % Verify
    % actual delay should correspond to desired delay
    actualDelay = finddelay(input,output);
    verifyEqual(testCase,actualDelay,delaySamples);
    
    % largest cross-correlation coefficient should be significantly larger
    % than second-largest
    sortedxCorrL = sort(xcorr(input(:,1), output(:,1)));
    xcorrRatioL = sortedxCorrL(end) / sortedxCorrL(end-1);
    verifyGreaterThan(testCase, xcorrRatioL, 9.5)

    sortedxCorrR = sort(xcorr(input(:,2), output(:,2)));
    xcorrRatioR = sortedxCorrR(end) / sortedxCorrR(end-1);
    verifyGreaterThan(testCase, xcorrRatioR, 9.5)
end

function testGammatoneFilterbankInAlgorithm(testCase)
% test if filterbank does near-perfect reconstruction of signal within the
% speech enhancement algorithm framework
    
    % Setup
    AlgorithmParameters = AlgorithmParametersConstructor();
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);
    % disable processing stage
    AlgorithmParameters.Cancellation = false;
    AlgorithmParameters.Enhancement = false;
    % extract expected delay
    delaySeconds = ...
        AlgorithmParameters.Gammatone.desiredDelayInSeconds;
    samplingrateHz = ...
        AlgorithmParameters.Gammatone.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    % create test signal
    input = rand(10000,2)-0.5;
    input = testInputSignal(input);
    
    % Exercise
    [output, ~] = ...
        speechEnhancement(input, AlgorithmParameters, AlgorithmStates);
    
    % Verify
    actualDelay = finddelay(input,output);
    for ch=1:2
        % actual delay should correspond to desired delay

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
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);
    % disable processing stage
    AlgorithmParameters.Cancellation = false;
    AlgorithmParameters.Enhancement = false;
    % extract expected delay
    delaySeconds = ...
        AlgorithmParameters.Gammatone.desiredDelayInSeconds;
    samplingrateHz = ...
        AlgorithmParameters.Gammatone.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    % create test signal
    input = rand(10000,2)-0.5;
    input = testInputSignal(input);

    % Exercise
    output = blockFeedingRoutine(input, BlockFeedingParameters,...
        AlgorithmParameters, AlgorithmStates);
    
    % Verify
    actualDelay = finddelay(input,output);
    for ch=1:2
        % actual delay should correspond to desired delay

        verifyEqual(testCase,actualDelay(ch),delaySamples);

        % largest cross-correlation coefficient should be significantly 
        % larger than second-largest
        sortedxCorr = sort(xcorr(input(:,ch),output(:,ch)));
        xcorrRatio = sortedxCorr(end)/sortedxCorr(end-1);
        verifyGreaterThan(testCase,xcorrRatio,9.5)
    end
end

function testGammatoneFilterbankInBlockFeederCustom(testCase)
% test if filterbank does near-perfect reconstruction of signal within the
% block feeding framework with custom block length
    
    % Setup
    BlockFeedingParameters.blockLength = 2000;
    AlgorithmParameters = AlgorithmParametersConstructor();
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);
    % disable processing stage
    AlgorithmParameters.Cancellation = false;
    AlgorithmParameters.Enhancement = false;
    % extract expected delay
    delaySeconds = ...
        AlgorithmParameters.Gammatone.desiredDelayInSeconds;
    samplingrateHz = ...
        AlgorithmParameters.Gammatone.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    % create test signal
    input = rand(10000,2)-0.5;
    input = testInputSignal(input);

    % Exercise
    output = blockFeedingRoutine(input,BlockFeedingParameters,...
        AlgorithmParameters, AlgorithmStates);
    
    % Verify
    actualDelay = finddelay(input,output);
    for ch=1:2
        % actual delay should correspond to desired delay

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
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);
    % disable processing stage
    AlgorithmParameters.Cancellation = false;
    AlgorithmParameters.Enhancement = false;
    % extract expected delay
    delaySeconds = ...
        AlgorithmParameters.Gammatone.desiredDelayInSeconds;
    samplingrateHz = ...
        AlgorithmParameters.Gammatone.samplingRateHz;
    delaySamples = round(delaySeconds * samplingrateHz);
    % create test signal
    input = rand(1000,2)-0.5;
    input = testInputSignal(input);

    % Exercise
    output = blockFeedingRoutine(input,BlockFeedingParameters,...
        AlgorithmParameters, AlgorithmStates);
    
    % Verify
    actualDelay = finddelay(input,output);
    for ch=1:2
        % actual delay should correspond to desired delay
        verifyEqual(testCase,actualDelay(ch),delaySamples);

        % largest cross-correlation coefficient should be significantly 
        % larger than second-largest
        sortedxCorr = sort(xcorr(input(:,ch),output(:,ch)));
        xcorrRatio = sortedxCorr(end)/sortedxCorr(end-1);
        verifyGreaterThan(testCase,xcorrRatio,7)
    end
end

function testGammatoneFilterbankInBlockFeederCompare(testCase)
% test if filterbank returns identical signals when run within the block
% feeding framework, running only one sample at a time, and when run as
% batch
    
    % Setup
    BlockFeedingParameters.blockLength = 1;
    AlgorithmParameters = AlgorithmParametersConstructor();
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);
    % disable processing stage
    AlgorithmParameters.Cancellation = false;
    AlgorithmParameters.Enhancement = false;

    AlgorithmStatesBlock = AlgorithmStates;
    AlgorithmStatesBatch = AlgorithmStates;
   
    % create test signal
    input = rand(1000,2)-0.5;
    input = testInputSignal(input);

    % Exercise
    blockOutput = blockFeedingRoutine(input, BlockFeedingParameters,...
        AlgorithmParameters, AlgorithmStatesBlock);
    [batchOutput, ~] = ...
        speechEnhancement(input, AlgorithmParameters, AlgorithmStatesBatch);
    
    % Verify
    verifyEqual(testCase,blockOutput,batchOutput,"AbsTol",1e-12)
end

% Periodicity Analysis
function testPeriodicityAnalysisExpectedSizes(testCase)
% test if the periodicityAnalysis function returns p0 candidate index
% values in the expected range and if the Sigma, Delta and SNR vectors all
% have the expected length

    % Setup
    AlgorithmParameters = AlgorithmParametersConstructor();
    
    % reduce number of frequency bands in the gammatone filterbank to
    % reduce computation time
    AlgorithmParameters.Gammatone.fHigh = 340;
    
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
        AlgorithmStatesConstructor(AlgorithmParameters);
    
    TestSignalParameters.nSamples = 1000;
    inputSignal = testSignalGenerator(TestSignalParameters);
    [subbandSignalArray, AlgorithmStates] = ...
        subbandDecompositionBinaural(inputSignal, AlgorithmStates);

    % Exercise
    [sigma, delta, snr, p0CandidateSampleIndexVectors, ~, ...
        p0SearchRangeSamplesVector] = periodicityAnalysis(...
        subbandSignalArray, AlgorithmParameters, AlgorithmStates);

    % Verify
    nP0Values = length(p0SearchRangeSamplesVector);

    nBands = AlgorithmParameters.Gammatone.nBands;
    for iBand=1:nBands
        verifyLessThanOrEqual(testCase,p0CandidateSampleIndexVectors.L{iBand},...
            nP0Values)
        verifyLessThanOrEqual(testCase,p0CandidateSampleIndexVectors.R{iBand},...
            nP0Values)
        
        verifyEqual(testCase,length(snr.L{iBand}),...
            length(find(p0CandidateSampleIndexVectors.L{iBand})))
        verifyEqual(testCase,length(sigma.L{iBand}),...
            length(find(p0CandidateSampleIndexVectors.L{iBand})))
        verifyEqual(testCase,length(delta.L{iBand}),...
            length(find(p0CandidateSampleIndexVectors.L{iBand})))
        verifyEqual(testCase,length(snr.R{iBand}),...
            length(find(p0CandidateSampleIndexVectors.R{iBand})))
        verifyEqual(testCase,length(sigma.R{iBand}),...
            length(find(p0CandidateSampleIndexVectors.R{iBand})))
        verifyEqual(testCase,length(delta.R{iBand}),...
            length(find(p0CandidateSampleIndexVectors.R{iBand})))
    end
end

function testCalcSigmaDeltaBinauralRangeSizesAndValues(testCase)
% test if Sigma/Delta range computation produces output of expected size
% (length-of-signal x length-of-p0-vector) with complex binaural signal as
% input, and if the returned FIFO equals the end of the input signal, and
% if the initial blocks of sigma and delta values (before the FIFO starts
% getting filled with past samples) are identical to the input signal since
% 0 is added/subtracted, and if the sigma and delta arrays equal the
% expected addition/subtraction of past samples from the current samples

    % Setup 
    testSignalReal = testSignalGenerator;
    testSignalImag = testSignalGenerator;
    subbandSignal.L = testSignalReal(:,1) + 1i*testSignalImag(:,1);
    subbandSignal.R = testSignalReal(:,2) + 1i*testSignalImag(:,2);

    p0Values = (32:160)';
    
    FilterStates.L.p0DetectionFIFO = zeros(1,max(p0Values));
    FilterStates.R.p0DetectionFIFO = zeros(1,max(p0Values));

    mode = 'range';
    
    % Exercise
    [sigma, delta, FilterStates] = calcSigmaDeltaBinaural(subbandSignal, ...
        p0Values, FilterStates, mode);
    
    % Verify
    expectedSize = [size(subbandSignal.L,1) length(p0Values)];
    verifySize(testCase,sigma.L,expectedSize)
    verifySize(testCase,sigma.R,expectedSize)
    verifySize(testCase,delta.L,expectedSize)
    verifySize(testCase,delta.R,expectedSize)

    verifyEqual(testCase,FilterStates.L.p0DetectionFIFO,...
        flip(subbandSignal.L(end-max(p0Values)+1:end).'))
    verifyEqual(testCase,FilterStates.R.p0DetectionFIFO,...
        flip(subbandSignal.R(end-max(p0Values)+1:end).'))

    for iPeriod = 1:length(p0Values)
        verifyEqual(testCase, sigma.L(1:min(p0Values), iPeriod), ...
            subbandSignal.L(1:min(p0Values)))
        verifyEqual(testCase, delta.L(1:min(p0Values), iPeriod), ...
            subbandSignal.L(1:min(p0Values)))
        verifyEqual(testCase, sigma.R(1:min(p0Values), iPeriod), ...
            subbandSignal.R(1:min(p0Values)))
        verifyEqual(testCase, delta.R(1:min(p0Values), iPeriod), ...
            subbandSignal.R(1:min(p0Values)))
    end
    for iSample = max(p0Values)+1:max(p0Values)+100
        for iPeriod = 1:length(p0Values)
            verifyEqual(testCase, sigma.L(iSample, iPeriod), ...
                subbandSignal.L(iSample) + ...
                subbandSignal.L(iSample-p0Values(iPeriod)))
            verifyEqual(testCase, delta.L(iSample, iPeriod), ...
                subbandSignal.L(iSample) - ...
                subbandSignal.L(iSample-p0Values(iPeriod)))
            verifyEqual(testCase, sigma.R(iSample, iPeriod), ...
                subbandSignal.R(iSample) + ...
                subbandSignal.R(iSample-p0Values(iPeriod)))
            verifyEqual(testCase, delta.R(iSample, iPeriod), ...
                subbandSignal.R(iSample) - ...
                subbandSignal.R(iSample-p0Values(iPeriod)))
        end
    end

end

function testCalcSigmaDeltaBinauralDiscreteSizes(testCase)
% test if Sigma/Delta discrete value computation produces output of 
% expected size (number of samples with detected p0) with complex 
% binaural signal as input, and if the returned FIFO equals the end of the 
% input signal, and if the computed Sigma and Delta values equal the
% expected addition/subtraction of past samples from the current samples

    % Setup 
    testSignalReal = testSignalGenerator;
    testSignalImag = testSignalGenerator;
    subbandSignal.L = testSignalReal(:,1) + 1i*testSignalImag(:,2);
    subbandSignal.R = testSignalReal(:,1) + 1i*testSignalImag(:,2);
    
    p0SearchRangeSamplesVector = (32:160)';

    FilterStates.L.p0CandidateFIFO = zeros(max(p0SearchRangeSamplesVector),1);
    FilterStates.R.p0CandidateFIFO = zeros(max(p0SearchRangeSamplesVector),1);
    
    p0CandidateSampleIndexVector.L = zeros(size(subbandSignal.L));
    p0CandidateSampleIndexVector.R = zeros(size(subbandSignal.R));

    p0CandidateSampleIndexVector.L(1200) = 107;
    p0CandidateSampleIndexVector.L(2000:2100) = 97;
    p0CandidateSampleIndexVector.R(1100) = 103;
    p0CandidateSampleIndexVector.R(2200:2400) = 82;

    mode = 'discrete';
    
    % Exercise
    [sigma, delta, FilterStates] = calcSigmaDeltaBinaural(subbandSignal, ...
        p0SearchRangeSamplesVector, FilterStates, mode, ...
        p0CandidateSampleIndexVector);
    
    % Verify
    expectedSize.L = [length(find(p0CandidateSampleIndexVector.L)) 1];
    expectedSize.R = [length(find(p0CandidateSampleIndexVector.R)) 1];
    verifySize(testCase,sigma.L,expectedSize.L)
    verifySize(testCase,sigma.R,expectedSize.R)
    verifySize(testCase,delta.L,expectedSize.L)
    verifySize(testCase,delta.R,expectedSize.R)

    verifyEqual(testCase,FilterStates.L.p0CandidateFIFO,...
        subbandSignal.L(end-max(p0SearchRangeSamplesVector)+1:end))
    verifyEqual(testCase,FilterStates.R.p0CandidateFIFO,...
        subbandSignal.R(end-max(p0SearchRangeSamplesVector)+1:end))

    verifyEqual(testCase, sigma.L(1), subbandSignal.L(1200) + ...
        subbandSignal.L(1200-p0SearchRangeSamplesVector(107)))
    verifyEqual(testCase, delta.L(1), subbandSignal.L(1200) - ...
        subbandSignal.L(1200-p0SearchRangeSamplesVector(107)))
    verifyEqual(testCase, sigma.L(2:102), subbandSignal.L(2000:2100) + ...
        subbandSignal.L((2000:2100)-p0SearchRangeSamplesVector(97)))
    verifyEqual(testCase, delta.L(2:102), subbandSignal.L(2000:2100) - ...
        subbandSignal.L((2000:2100)-p0SearchRangeSamplesVector(97)))
    verifyEqual(testCase, sigma.R(1), subbandSignal.R(1100) + ...
        subbandSignal.R(1100-p0SearchRangeSamplesVector(103)))
    verifyEqual(testCase, delta.R(1), subbandSignal.R(1100) - ...
        subbandSignal.R(1100-p0SearchRangeSamplesVector(103)))
    verifyEqual(testCase, sigma.R(2:202), subbandSignal.R(2200:2400) + ...
        subbandSignal.R((2200:2400)-p0SearchRangeSamplesVector(82)))
    verifyEqual(testCase, delta.R(2:202), subbandSignal.R(2200:2400) - ...
        subbandSignal.R((2200:2400)-p0SearchRangeSamplesVector(82)))
end

function testAbsoluteSquareLPFilterBinauralSizes(testCase)
% test if the absoluteSquareLPFilterBinaural function generates output
% structs of equal size to the input and if the returned filter states are
% non-zero, and if the output values are >= 0 as they should be due to
% absolute squaring

    % Setup
    sigmaIn.L = testSignalGenerator;
    sigmaIn.R = testSignalGenerator;
    deltaIn.L = testSignalGenerator;
    deltaIn.R = testSignalGenerator;
    AlgorithmParameters = AlgorithmParametersConstructor();
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
        AlgorithmStatesConstructor(AlgorithmParameters);
    FilterStates.L = AlgorithmStates.L.ProcessingStates{1};
    FilterStates.R = AlgorithmStates.R.ProcessingStates{1};
    % Exercise
    [sigmaOut, deltaOut, FilterStates] = ...
        absoluteSquareLPFilterBinaural(sigmaIn, deltaIn, AlgorithmParameters, FilterStates);

    % Verify
    verifyNotEqual(testCase,FilterStates.L.sigmaNormLP,0);
    verifyNotEqual(testCase,FilterStates.R.sigmaNormLP,0);
    verifyNotEqual(testCase,FilterStates.L.deltaNormLP,0);
    verifyNotEqual(testCase,FilterStates.R.deltaNormLP,0);
    expectedSize = size(sigmaIn.L);
    verifySize(testCase,deltaIn.R,expectedSize)
    verifySize(testCase,sigmaOut.L,expectedSize)
    verifySize(testCase,sigmaOut.R,expectedSize)
    verifySize(testCase,deltaOut.L,expectedSize)
    verifySize(testCase,deltaOut.R,expectedSize)
    verifyGreaterThanOrEqual(testCase, sigmaOut.L, 0)
    verifyGreaterThanOrEqual(testCase, sigmaOut.R, 0)
    verifyGreaterThanOrEqual(testCase, deltaOut.L, 0)
    verifyGreaterThanOrEqual(testCase, deltaOut.R, 0)
end

function testFirstOrderLowPassSizes(testCase)
% test if the first order lowpass filter function returns output signal of
% equal size to input signal and returns correct filter state for a
% single-channel input signal

    % Setup
    AlgorithmParameters = AlgorithmParametersConstructor();
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
        AlgorithmStatesConstructor(AlgorithmParameters);
    filterState = AlgorithmStates.L.ProcessingStates{1}.sigmaNormLP;
    inputSignal = testSignalGenerator;
    inputSignal = inputSignal(:,1);
    
    % Exercise
    [outputSignal, filterState] = firstOrderLowPass(inputSignal, ...
        filterState, AlgorithmParameters);

    % Verify
    expectedSize = size(inputSignal);
    verifySize(testCase,outputSignal,expectedSize)
    verifyEqual(testCase,filterState,outputSignal(end,:))
end

function testFirstOrderLowPassBinauralSizes(testCase)
% test if the first order lowpass filter function returns output signal of
% equal size to input signal and returns correct filter state for a
% two-channel input signal

    % Setup
    AlgorithmParameters = AlgorithmParametersConstructor();
    filterState = [0 0];
    inputSignal = testSignalGenerator;
    
    % Exercise
    [outputSignal, filterState] = firstOrderLowPass(inputSignal, ...
        filterState, AlgorithmParameters);

    % Verify
    expectedSize = size(inputSignal);
    verifySize(testCase,outputSignal,expectedSize)
    verifyEqual(testCase,filterState,outputSignal(end,:))
end

function testFirstOrderLowPassComplexSizes(testCase)
% test if the first order lowpass filter function returns output signal of
% equal size to input signal and returns correct filter state for a
% single-channel complex input signal

    % Setup
    AlgorithmParameters = AlgorithmParametersConstructor();
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
        AlgorithmStatesConstructor(AlgorithmParameters);
    filterState = AlgorithmStates.L.ProcessingStates{1}.sigmaNormLP;
    testSignal = testSignalGenerator;
    inputSignal = testSignal(:,1) + 1i*testSignal(:,2);

    % Exercise
    [outputSignal, filterState] = firstOrderLowPass(inputSignal, ...
        filterState, AlgorithmParameters);

    % Verify
    expectedSize = size(inputSignal);
    verifySize(testCase,outputSignal,expectedSize)
    verifyEqual(testCase,filterState,outputSignal(end,:))
end
    
function testFirstOrderLowpassFilterStateRelay(testCase)
% test if block-wise filtering with first order lowpass results in the same
% output as batch filtering the entire signal, thus testing the function of
% the filter state variable
    % Setup
    AlgorithmParameters = AlgorithmParametersConstructor();
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
        AlgorithmStatesConstructor(AlgorithmParameters);
    filterState = AlgorithmStates.L.ProcessingStates{1}.sigmaNormLP;
    inputSignal = testSignalGenerator;
    inputSignal = inputSignal(:,1);

    blockLength = 1000;
    nBlocks = floor(size(inputSignal,1)/blockLength);
    
    % Exercise
    [batchOutputSignal, ~] = firstOrderLowPass(inputSignal, ...
        filterState, AlgorithmParameters);
    
    blockOutputSignal = zeros(size(inputSignal));
    for iBlock=1:nBlocks
        signalIndices = round((1:blockLength)+(iBlock-1)*blockLength);
        [blockOutputSignal(signalIndices), filterState] = ...
            firstOrderLowPass(inputSignal(signalIndices), filterState, ...
            AlgorithmParameters);
    end

    % Verify
    verifyEqual(testCase,blockOutputSignal,batchOutputSignal)
end

function testCalcSNRBinauralSizes(testCase)
% test if a binaural set of input signals in the calcSNRBinaural funtion
% generates output structures of the same size

    % Setup 
    sigma.L = testSignalGenerator;
    sigma.R = testSignalGenerator;
    delta.L = testSignalGenerator;
    delta.R = testSignalGenerator;

    % Exercise
    snr = calcSNRBinaural(sigma, delta);

    % Verify
    expectedSize = size(sigma.L);
    verifySize(testCase,snr.L,expectedSize)
    verifySize(testCase,snr.R,expectedSize)
end

function testCalcSNRBinauralExpectedValue(testCase)
% test if a binaural set of input signals in the calcSNRBinaural funtion
% generates output of correct value

    % Setup 
    sigmaVal = 3;
    deltaVal = 4;

    sigma.L = sigmaVal*ones(100,1);
    sigma.R = sigmaVal*ones(100,1);
    delta.L = deltaVal*ones(100,1);
    delta.R = deltaVal*ones(100,1);

    % Exercise
    snr = calcSNRBinaural(sigma, delta);

    % Verify
    expectedValue = (sigmaVal/deltaVal)*ones(100,1);
    verifyEqual(testCase, snr.L, expectedValue, "AbsTol", 1e-15)
    verifyEqual(testCase, snr.R, expectedValue, "AbsTol",1e-15)
end

function testSubbandSnrPeakDetectionBinauralSizes(testCase)
% test if the subbandSnrPeakDetectionBinaural function returns arrays of
% the expected dimension and if they only contain values less or equal to
% the maximum p0 value in samples

    % Setup
    nP0Values = 50;
    for iP0Value = 1:nP0Values
        testSignal = testSignalGenerator;
        snr.L(:,iP0Value) = testSignal(:,1);
        testSignal = testSignalGenerator;
        snr.R(:,iP0Value) = testSignal(:,2);
    end
    % Exercise
    p0Candidates = subbandSnrPeakDetectionBinaural(snr);

    % Validate
    expectedSize = [size(snr.L,1), 1];
    verifySize(testCase, p0Candidates.L, expectedSize)
    verifySize(testCase, p0Candidates.R, expectedSize)

    verifyLessThanOrEqual(testCase, p0Candidates.L, nP0Values)
    verifyLessThanOrEqual(testCase, p0Candidates.R, nP0Values)
end

function testSubbandSnrPeakDetectionBinauralValues(testCase)
% test if the subbandSnrPeakDetectionBinaural function returns the expected
% values, i.e. the indices of the maximum values along the correct axis

    % Setup
    nSamples = 30;
    for iSample = 1:nSamples
        snr.L(iSample,:) = sin(pi/100:pi/100:10*pi);
        snr.R(iSample,:) = (linspace(0,1,1000)).*sin(pi/100:pi/100:10*pi);
    end

    % Exercise
    p0Candidates = subbandSnrPeakDetectionBinaural(snr);

    % Validate
%     expectedP0Candidates.L = 50*ones(nSamples,1);
%     expectedP0Candidates.L = 851*ones(nSamples,1);
    [~, maxIndices.L] = max(snr.L, [], 2);
    [~, maxIndices.R] = max(snr.R, [], 2);

    verifyEqual(testCase, p0Candidates.L, maxIndices.L)
    verifyEqual(testCase, p0Candidates.R, maxIndices.R)
end

% Direction-of-arrival estimation
function testIldBinauralLowPassFilter1000Hz(testCase)
% test if the function effectively filters out a frequency of 1000 Hz

    % Setup
    samplingRateHz = 10000;
    dt = 1/samplingRateHz;
    
    timeVectorSeconds = dt:dt:1;
    freqHz = 1000;
    signal = sin(2*pi*freqHz*timeVectorSeconds);
    subbandSignal.L = signal;
    subbandSignal.R = signal;

    state = [0; 0];
    States.L.ildLP = state;
    States.R.ildLP = state;

    % Exercise
    [subbandSignalLp, ~] = ...
        ildBinauralLowPassFilter(States, subbandSignal, samplingRateHz);

    % Verify
    verifyLessThan(testCase, max(abs(subbandSignalLp.L))/...
        max(abs(subbandSignal.L)), 0.02)
    verifyLessThan(testCase, max(abs(subbandSignalLp.R))/...
        max(abs(subbandSignal.R)), 0.02)
end

function testIldBinauralLowPassFilter20Hz(testCase)
% test if the function leaves a frequency of 10 Hz untouched

    % Setup
    samplingRateHz = 10000;
    dt = 1/samplingRateHz;
    
    timeVectorSeconds = dt:dt:1;
    freqHz = 10;
    signal = sin(2*pi*freqHz*timeVectorSeconds);
    subbandSignal.L = signal;
    subbandSignal.R = signal;

    state = [0; 0];
    States.L.ildLP = state;
    States.R.ildLP = state;

    % Exercise
    [subbandSignalLp, ~] = ...
        ildBinauralLowPassFilter(States, subbandSignal, samplingRateHz);

    % Verify
    verifyGreaterThan(testCase, max(abs(subbandSignalLp.L))/...
        max(abs(subbandSignal.L)), 0.98)
    verifyGreaterThan(testCase, max(abs(subbandSignalLp.R))/...
        max(abs(subbandSignal.R)), 0.98)
end

function testIldBinauralLowPassFilterStateRelay(testCase)
% test if block-wise filtering with lowpass for ILD results in the same
% output as batch filtering the entire signal, thus testing the function of
% the filter state variable
    % Setup
    samplingRateHz = 16000;
    state = [0; 0];
    BatchStates.L.ildLP = state;
    BatchStates.R.ildLP = state;
    BlockStates.L.ildLP = state;
    BlockStates.R.ildLP = state;
    
    testSignal = testSignalGenerator;
    inputSignal.L = testSignal(:,1);
    inputSignal.R = testSignal(:,2);

    blockLength = 1000;
    nBlocks = floor(size(inputSignal.L,1)/blockLength);
    
    % Exercise
    [batchOutputSignalLp, BatchStates] = ...
        ildBinauralLowPassFilter(BatchStates, inputSignal, samplingRateHz);
    
    blockOutputSignalLp.L = zeros(size(inputSignal.L));
    blockOutputSignalLp.R = zeros(size(inputSignal.R));
    for iBlock=1:nBlocks

        signalIndices = round((1:blockLength)+(iBlock-1)*blockLength);

        blockIn.L = inputSignal.L(signalIndices);
        blockIn.R = inputSignal.R(signalIndices);

        [blockOut, BlockStates] = ...
            ildBinauralLowPassFilter(BlockStates, blockIn, samplingRateHz);

        blockOutputSignalLp.L(signalIndices) = blockOut.L;
        blockOutputSignalLp.R(signalIndices) = blockOut.R;

    end

    % Verify
    verifyEqual(testCase, blockOutputSignalLp, batchOutputSignalLp, ...
        "AbsTol", 1e-14)
    verifyEqual(testCase, BlockStates, BatchStates, "AbsTol", 1e-14)
end

function testDisambiguateIpdActive(testCase)
% test if function sets the sign as expected, when abs(ILD) is greater than
% 2.5 dB

    % Setup
    ildVec = linspace(0,2*pi);
    ildDb = 5*square(ildVec);

    inIpdRad = -1*square(ildVec);

    % Exercise
    outIpdRad = disambiguateIpd(inIpdRad, ildDb);

    % Verify
    verifyEqual(testCase, -inIpdRad, outIpdRad)

end

function testDisambiguateIpdInactive(testCase)
% test if function leaves the sign as expected, when abs(ILD) is less than
% 2.5 dB

    % Setup
    ildVec = linspace(0,2*pi);
    ildDb = 1*square(ildVec);

    inIpdRad = -1*square(ildVec);

    % Exercise
    outIpdRad = disambiguateIpd(inIpdRad, ildDb);

    % Verify
    verifyEqual(testCase, inIpdRad, outIpdRad)

end

function testCalcIvs(testCase)

% test if the IVS mask calculation returns a sufficiently low count of
% samples passed for incoherent signals and sufficiently high for coherent
% signals

% % non-coherent noise
% binauralSignal(:,1) = randn(10000,1) + 1i*randn(10000,1);
% binauralSignal(:,2) = randn(10000,1) + 1i*randn(10000,1);
% verifyLessThan(testCase, ratio, 0.0005)
% % coherent noise 0°
% signal = randn(10000,1) + 1i*randn(10000,1);
% hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
% binauralSignal = SOFAspat(signal, hrtf, 0, 0);
% verifyGreaterThan(testCase, ratio, 0.82)
% % coherent noise -60°
% signal = randn(10000,1) + 1i*randn(10000,1);
% hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
% binauralSignal = SOFAspat(signal, hrtf, -60, 0);
% verifyGreaterThan(testCase, ratio, 0.82)
% 
% % non-coherent speech
% signal = audioread('p298_097.wav');
% binauralSignal(:,1) = SOFAspat(signal, hrtf, 40, 0);
% binauralSignal(:,2) = SOFAspat(signal, hrtf, -60, 0);
% verifyLessThan(testCase, ratio, 0.0005)
% % coherent speech 0°
% signal = audioread('p298_097.wav');
% hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
% binauralSignal = SOFAspat(signal, hrtf, 0, 0);
% verifyGreaterThan(testCase, ratio, 0.0005)
% % coherent speech -60°
% signal = audioread('p298_097.wav');
% hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
% binauralSignal = SOFAspat(signal, hrtf, -60, 0);
% verifyGreaterThan(testCase, ratio, 0.0005)
% 
%     % Setup
%     ratio=zeros(10000,1);
%     for i=1:10000
%     signal = randn(10000,1)+1i*randn(10000,1);
% %     signal = audioread('p298_097.wav');
% %     hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
%     binauralSignal = SOFAspat(signal, hrtf, -60, 0);
%     itf = binauralSignal(:,1) .* conj(binauralSignal(:,2));
% 
% %     States = AlgorithmStatesConstructor(AlgorithmParametersConstructor());
% %     States = States.Binaural.ProcessingStates{1};
% 
%     % Exercise
%     ivs = calcIvs(itf, 0.04, 44100, 0.98, States);
% 
%     % Verify
%     ratio(i) = numel(nonzeros(ivs))/numel(ivs);
%     end
%     verifyGreaterThan(testCase, ratio, 0.6)
end

function testCalcIvsFilterStateRelay(testCase)
% test if block-wise calculation of the IVS mask results output as batch 
% processing the entire signal, thus testing the function of
% the filter states and previous value in the States struct
    % Setup
    signal = audioread('p298_097.wav');
    hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
    binauralSignal = SOFAspat(signal, hrtf, -60, 0);
    binauralSignal = binauralSignal(1:10000,:);
    itf = binauralSignal(:,1) .* conj(binauralSignal(:,2));

    States = AlgorithmStatesConstructor(AlgorithmParametersConstructor());
    States = States.Binaural.ProcessingStates{1};

    % Exercise
    ivsMaskBatch = calcIvs(itf, 0.04, 10000, 0.98, States);
    
    for iBlock = 1:100
        [ivsMaskBlock((1:100)+(iBlock-1)*100,1), States] = ...
            calcIvs(itf((1:100)+(iBlock-1)*100), 0.04, 10000, 0.98, States);
    end

    % Verify
    verifyEqual(testCase, ivsMaskBlock, ivsMaskBatch)
end

% Azimuth estimation
function testAzimuthEstimation(testCase)
end

function testInterauralToAzimuthMapping(testCase)
end

% Harmonic enhancement
function testHarmonicEnhancemne(testCase)
%     [subbandSignalArray.L, targetSampleIndices.L{iBand}, ...
%         interfererSampleIndices.L{iBand}] = harmonicEnhancement(...
%         subbandSignalArray.L, iBand, ivsMask, ...
%         p0DetectedIndexVectors.L, azimuthDegCells, snrDesired.L, ...
%         sigmaDesired.L, deltaDesired.L, AlgorithmParameters);
end
%% Full Speech enhancement algorithm
function testSpeechEnhancementInBlockFeederCompare(testCase)
% test if speech enhancement returns identical output signal and filter
% states when run within the block feeding framework, running only one 
% sample at a time, running one block at a time, and when run as batch
    
    % Setup
    AlgorithmParameters = AlgorithmParametersConstructor();
    
%     AlgorithmParameters.lookuptableType = 'ipd';
    
    hrtf = SOFAload('HRIR_KEMAR_DV0001_3.sofa',[5 2],'R');
    AlgorithmParameters.Gammatone.samplingRateHz = hrtf.Data.SamplingRate;
    AlgorithmParameters.lookuptable = ...
        interauralToAzimuthLookuptable(hrtf, AlgorithmParameters);

    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
        AlgorithmStatesConstructor(AlgorithmParameters);

    SampleFeedingParameters.blockLength = 1;
    BlockFeedingParameters.blockLength = 100;

    SampleStates = AlgorithmStates;
    BlockStates = AlgorithmStates;
    BatchStates = AlgorithmStates;
   
    % create test signal
    input = rand(500,2)-0.5;
    input = testInputSignal(input);

    % Exercise
    fprintf('\n sample-by-sample processing\n')
    tic
    [sampleOutput, SampleStates] = blockFeedingRoutine(input,...
        SampleFeedingParameters, AlgorithmParameters, SampleStates);
    toc
    fprintf('\n block-by-block processing\n')
    tic
    [blockOutput, BlockStates] = blockFeedingRoutine(input,...
        BlockFeedingParameters, AlgorithmParameters, BlockStates);
    toc
    fprintf('\n batch processing\n')
    tic
    [batchOutput, BatchStates] = ...
        speechEnhancement(input, AlgorithmParameters, BatchStates);
    toc

    % Verify
    verifyEqual(testCase,SampleStates,BatchStates,"AbsTol",1e-12)
    verifyEqual(testCase,BlockStates,BatchStates,"AbsTol",1e-12)
    verifyEqual(testCase,sampleOutput,batchOutput,"AbsTol",1e-14)
    verifyEqual(testCase,blockOutput,batchOutput,"AbsTol",1e-14)
    

    
end

%% LP filter testing
function testFirstOrderLPFilterAgainstFilter(testCase)
% verify that the firstOrderLPFilter.m function functions equally to the
% MATLAB-internal filter function

% Setup
AlgorithmParameters = AlgorithmParametersConstructor;
fs = AlgorithmParameters.Gammatone.samplingRateHz;
T = 1/fs;
tau = 0.1;
a = exp(-(T/tau));

inputSignal = testSignalGenerator;

% Exercise
[ownOutput, ~] = firstOrderLowPass(inputSignal, 0, AlgorithmParameters, tau);

matlabOutput = filter([1-a], [1, -a], inputSignal);

% Verify
verifyEqual(testCase, ownOutput, matlabOutput)
end