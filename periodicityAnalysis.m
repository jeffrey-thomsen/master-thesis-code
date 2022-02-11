% Periodicity analysis stage of the Thomsen2022 speech enhancement
% algorithm. Uses a set of binaural signals that has been decomposed by a
% gammatone filterbank and returns Sigma, Delta and SNR values for each
% subband sample at which a periodic component was detected. The frequency
% range for the periodicity analysis is defined in p0SearchRange
function [sigmaCells, deltaCells, snrCells, ...
  p0CandidateSampleIndexVectorCells, AlgorithmStates] = ...
  periodicityAnalysis(subbandSignalArray, AlgorithmParameters, AlgorithmStates)

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
        FilterStates.L = AlgorithmStates.L.ProcessingStates{iBand};
        FilterStates.R = AlgorithmStates.R.ProcessingStates{iBand};
        
        % calculate normalized SNR
        [normalizedSigma, normalizedDelta] = ...
            calcSigmaDeltaBinaural(normalizedSubbandSignal, ...
                p0SearchRangeSamplesVector, 'range');
    
        [normalizedSigma, normalizedDelta, FilterStates] = ...
            absoluteSquareLPFilterBinaural(normalizedSigma, ...
            normalizedDelta, AlgorithmParameters, FilterStates);
    
        normalizedSnr = calcSNRBinaural(normalizedSigma, normalizedDelta);
        
        % intra-subband SNR peak detection
        p0CandidateSampleIndexVector = ...
            subbandSnrPeakDetectionBinaural(normalizedSnr);

        % calculate Sigmas, Deltas and SNR for p0 candidates
        [p0CandidateSigma, p0CandidateDelta] = ...
            calcSigmaDeltaBinaural(subbandSignal, ...
            p0SearchRangeSamplesVector, 'discrete', ...
            p0CandidateSampleIndexVector);
    
        p0CandidateSnr = calcSNRBinaural(p0CandidateSigma, p0CandidateDelta);

        % write into structs
        AlgorithmStates.L.AlgorithmStates{iBand} = FilterStates.L;
        AlgorithmStates.R.AlgorithmStates{iBand} = FilterStates.R;
        
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