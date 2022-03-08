% Periodicity analysis stage of the Thomsen2022 speech enhancement
% algorithm. Uses a set of binaural signals that has been decomposed by a
% gammatone filterbank and returns Sigma, Delta and SNR values for each
% subband sample at which a periodic component was detected. The frequency
% range for the periodicity analysis is defined in p0SearchRange
%
% N - length of signal
% M - number of subbands
%
% Input:
% subbandSignalArray - struct containing NxM arrays of complex-valued
% subband samples from a gammatone filterbank, for left and right channel
% AlgorithmParameters - struct containing all simulation parameters
% AlgorithmStates - struct containing filter states and FIFO arrays that
% need to be handed over to the next signal block/sample to be processed
%
% Output:
% sigmaCells, deltaCells, snrCells - cells of size M containing Sigma, 
% Delta and SNR values for each detected periodic sample in each subband
% p0CandidateSampleIndexVectorCells - cells of size M specifying at the
% positions of detected periodic samples within each subband signal block
% the index pointing to the detected period
% AlgorithmStates - see above
% p0SearchRangeSamplesVector - range vector for mapping the indices in
% p0CandidateSampleIndexVectorCells to the actual number of samples that
% the detected period corresponds to
function [sigmaCells, deltaCells, snrCells, ...
  p0CandidateSampleIndexVectorCells, AlgorithmStates, ...
  p0SearchRangeSamplesVector] = periodicityAnalysis(subbandSignalArray, ...
  AlgorithmParameters, AlgorithmStates)

    % convert p0 frequency search range in Hz into number of samples
    p0SearchRangeHz = AlgorithmParameters.p0SearchRangeHz;
    samplingRateHz = AlgorithmParameters.Gammatone.samplingRateHz;

    nMinSamplesP0Detection = floor(samplingRateHz/p0SearchRangeHz(2));
    nMaxSamplesP0Detection = ceil(samplingRateHz/p0SearchRangeHz(1));
    p0SearchRangeSamplesVector = ...
        (nMinSamplesP0Detection:nMaxSamplesP0Detection)';

    % normalize subband signals
    normalizedSubbandSignalArray.L = sign(subbandSignalArray.L);
    normalizedSubbandSignalArray.R = sign(subbandSignalArray.R);

    % preallocate data structures for results
    nBands = AlgorithmParameters.Gammatone.nBands;
    p0CandidateSampleIndexVectorCells.L = cell(1,nBands);
    p0CandidateSampleIndexVectorCells.R = cell(1,nBands);
    sigmaCells.L = cell(1,nBands);
    sigmaCells.R = cell(1,nBands);
    deltaCells.L = cell(1,nBands);
    deltaCells.R = cell(1,nBands);
    snrCells.L = cell(1,nBands);
    snrCells.R = cell(1,nBands);

    for iBand = 1:nBands
        subbandSignal.L = subbandSignalArray.L(:,iBand);
        subbandSignal.R = subbandSignalArray.R(:,iBand);
        normalizedSubbandSignal.L = normalizedSubbandSignalArray.L(:,iBand);
        normalizedSubbandSignal.R = normalizedSubbandSignalArray.R(:,iBand);
        States.L = AlgorithmStates.L.ProcessingStates{iBand};
        States.R = AlgorithmStates.R.ProcessingStates{iBand};

        if AlgorithmParameters.ChenP0Detection
            % like in Chen2015, use envelope for center frequencies above
            % 1.5kHz, use regular output below.
            % This is an alternative to the method proposed by Bruemann.
            if AlgorithmStates.L.GammatoneStates.analyzer.center_frequencies_hz(iBand)>1500
                normalizedSubbandSignal.L = abs(subbandSignal.L);
                normalizedSubbandSignal.R = abs(subbandSignal.R);
            else
                normalizedSubbandSignal.L = subbandSignal.L;
                normalizedSubbandSignal.R = subbandSignal.R;
            end
        end

        % compute sigma and delta values for all p0 in search range from
        % normalized subband signal
        [normalizedSigma, normalizedDelta, States] = ...
            calcSigmaDeltaBinaural(normalizedSubbandSignal, ...
            p0SearchRangeSamplesVector, States, 'range');
        
        % pre-process normalized sigmas and deltas for SNR computation
        [normalizedSigma, normalizedDelta, States] = ...
            absoluteSquareLPFilterBinaural(normalizedSigma, ...
            normalizedDelta, AlgorithmParameters, States);
        
        % compute normalized SNR for every p0 in search range
        normalizedSnr = calcSNRBinaural(normalizedSigma, normalizedDelta);
        
        % intra-subband SNR peak detection: find p0 candidate for each
        % signal sample
        p0CandidateSampleIndexVector = ...
            subbandSnrPeakDetectionBinaural(normalizedSnr);

        % calculate Sigmas, Deltas and SNR for p0 candidates from subband
        % signal
        [p0CandidateSigma, p0CandidateDelta, States] = ...
            calcSigmaDeltaBinaural(subbandSignal, ...
            p0SearchRangeSamplesVector, States, 'discrete', ...
            p0CandidateSampleIndexVector);
        p0CandidateSnr = calcSNRBinaural(p0CandidateSigma, p0CandidateDelta);

        % write states back into global structs
        AlgorithmStates.L.ProcessingStates{iBand} = States.L;
        AlgorithmStates.R.ProcessingStates{iBand} = States.R;
        
        % write computed values into structs for enhancement stage
        p0CandidateSampleIndexVectorCells.L{iBand} = ...
            p0CandidateSampleIndexVector.L;
        p0CandidateSampleIndexVectorCells.R{iBand} = ...
            p0CandidateSampleIndexVector.R;
        
        sigmaCells.L{iBand} = p0CandidateSigma.L;
        sigmaCells.R{iBand} = p0CandidateSigma.R;
        deltaCells.L{iBand} = p0CandidateDelta.L;
        deltaCells.R{iBand} = p0CandidateDelta.R;
        snrCells.L{iBand} = p0CandidateSnr.L;
        snrCells.R{iBand} = p0CandidateSnr.R;

    end

end