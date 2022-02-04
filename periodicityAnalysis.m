function [sigma, delta, snr] = ...
    periodicityAnalysis(subbandSignals, AlgorithmParameters)

    % convert p0 frequency search range to number of samples vector
    p0searchRangeHz = AlgorithmParameters.p0SearchRangeHz;
    samplingRateHz = AlgorithmParameters.GammatoneParameters.samplingRateHz;

    nMinSamplesP0Detect = floor(samplingRateHz/p0searchRangeHz(2));
    nMaxSamplesP0Detect = ceil(samplingRateHz/p0searchRangeHz(1));
    p0DetectSamplesVector = nMinSamplesP0Detect:nMaxSamplesP0Detect;
    nP0DetectSamples = length(p0DetectSamplesVector);

    % calculate normalized SNR
    norm.L = sign(subbandSignals.L);
    norm.R = sign(subbandSignals.R);
    
    [norm.sigma, norm.delta] = ...
        calcSigmaDelta(norm, p0DetectSamplesVector);

    norm = absoluteSquare(norm);

    norm = LPFilter(norm);

    norm.SNR = calcSNR(norm);
    
    % intra-subband SNR peak detection
    p0candidates = subbandSnrPeakDetection(norm.SNR);

    % calculate Sigmas, Deltas and SNR for p0 candidates
    [p0Detected.sigma, p0Detected.delta] = ...
        calcSigmaDelta(subbandSignals, p0candidates);

    p0Detected.SNR = calcSNR(p0Detected);



    snr = subbandSignals;
    sigma = snr;
    delta = sigma;
end