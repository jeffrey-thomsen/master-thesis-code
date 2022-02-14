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
    [sigma, delta, snr, p0CandidateSampleIndexVectors, ~] =...
        periodicityAnalysis(subbandSignalArray, AlgorithmParameters, ...
        AlgorithmStates);

    % Verify
    p0SearchRangeHz = AlgorithmParameters.p0SearchRangeHz;
    samplingRateHz = AlgorithmParameters.Gammatone.samplingRateHz;
    nMinSamplesP0Detection = floor(samplingRateHz/p0SearchRangeHz(2));
    nMaxSamplesP0Detection = ceil(samplingRateHz/p0SearchRangeHz(1));
    p0SearchRangeSamplesVector = ...
        (nMinSamplesP0Detection:nMaxSamplesP0Detection)';
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

function testPeriodicityAnalysisInBlockFeederCompare(testCase)
% test if periodicity analysis returns identical output signal and filter
% states when run within the block feeding framework, running only one 
% sample at a time, and when run as batch
    
    % Setup
    BlockFeedingParameters.blockLength = 1;
    AlgorithmParameters = AlgorithmParametersConstructor();
    [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);

    BlockStates = AlgorithmStates;
    BatchStates = AlgorithmStates;
   
    % create test signal
    input = rand(500,2)-0.5;
    input = testInputSignal(input);

    % Exercise
    [blockOutput, BlockStates] = blockFeedingRoutine(input,...
        BlockFeedingParameters, AlgorithmParameters, BlockStates);
    [batchOutput, BatchStates] = ...
        speechEnhancement(input, AlgorithmParameters, BatchStates);

    % Verify
    verifyEqual(testCase,blockOutput,batchOutput,"AbsTol",1e-12)
    
    verifyEqual(testCase,BlockStates,BatchStates,"AbsTol",1e-12)
    
end

function testCalcSigmaDeltaBinauralRangeSizes(testCase)
% test if Sigma/Delta range computation produces output of expected size
% (length-of-signal x length-of-p0-vector) with complex binaural signal as
% input, and if the returned FIFO equals the end of the input signal

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
end

function testCalcSigmaDeltaBinauralDiscreteSizes(testCase)
% test if Sigma/Delta discrete value computation produces output of 
% expected size (number of samples with detected p0) with complex 
% binaural signal as input, and if the returned FIFO equals the end of the 
% input signal

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
end

function testAbsoluteSquareLPFilterBinauralSizes(testCase)
% test if the absoluteSquareLPFilterBinaural function generates output
% structs of equal size to the input and if the returned filter states are
% non-zero

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
% test if a binaural set of input signals in the calsSNRBinaural funtion
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

function testSubbandSnrPeakDetectionBinauralSizes(testCase)
% test if the subbandSnrPeakDetectionBinaural function returns arrays of
% the expected dimension and if they only contain values less or equal to
% the maximum p0 value in samples

    % Setup
    nP0Values = 50;
    for iP0Value = 1:nP0Values
        testSignal = testSignalGenerator;
        snr.L(:,iP0Value) = testSignal(:,1) + 1i*testSignal(:,2);
        testSignal = testSignalGenerator;
        snr.R(:,iP0Value) = testSignal(:,1) + 1i*testSignal(:,2);
    end
    % Exercise
    p0Candidates = subbandSnrPeakDetectionBinaural(snr);

    % Validate
    expectedSize = [size(snr.L,1), 1];
    verifySize(testCase,p0Candidates.L,expectedSize)
    verifySize(testCase,p0Candidates.R,expectedSize)

    verifyLessThanOrEqual(testCase,p0Candidates.L,nP0Values)
    verifyLessThanOrEqual(testCase,p0Candidates.R,nP0Values)
end

% Direction-of-arrival estimation