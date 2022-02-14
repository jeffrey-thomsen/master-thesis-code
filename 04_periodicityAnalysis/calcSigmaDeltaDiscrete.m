function [sigma, delta, paddingFIFO] = calcSigmaDeltaDiscrete(subbandSignal, ...
  p0SearchRangeSamplesVector, p0CandidateIndexVector, paddingFIFO)
    % Computes sigmas and deltas for specific samples of the subbandSignal,
    % each at a specific value in the past.
    % Specified in p0CandidateIndexVector are the indices pointing at that
    % value in the past within the p0SearchRangeSamplesVector, stored at
    % the position corresponding to the signal sample index of
    % subbandSignal for which the specific sigma and delta are computed.

    sigmaFunc=@(a,b) bsxfun(@plus, a,b);
    deltaFunc=@(a,b) bsxfun(@minus,a,b);
    pastFunc=@(s,p0Indices,p0Values,shift) ...
        s(bsxfun(@minus,p0Indices+shift,p0Values));
    
    paddingFIFOshift = length(paddingFIFO);

    %% at each p0 candidate index compute sigma+delta with the p0CandidateSample
    % indices of signal samples at which a p0 candidate was detected
    p0CandidateSignalSampleIndexVector = find(p0CandidateIndexVector);
    % indices of corresponding p0 candidate values 
    p0CandidateNonZeroIndexVector = ...
        p0CandidateIndexVector(p0CandidateSignalSampleIndexVector);
    % map indices to p0 candidate values (in samples)
    p0CandidateValueVector = p0SearchRangeSamplesVector(p0CandidateNonZeroIndexVector);
    
    PaddedSubbandSignal = [paddingFIFO; subbandSignal];
    pastSamplesAtP0Candidate = pastFunc(PaddedSubbandSignal, ...
        p0CandidateSignalSampleIndexVector, p0CandidateValueVector, paddingFIFOshift);

    sigma = sigmaFunc(subbandSignal(p0CandidateSignalSampleIndexVector), ...
        pastSamplesAtP0Candidate);
    delta = deltaFunc(subbandSignal(p0CandidateSignalSampleIndexVector), ...
        pastSamplesAtP0Candidate);

    %% update padding FIFO
    if paddingFIFOshift<length(subbandSignal)
        paddingFIFO = subbandSignal(end-paddingFIFOshift+1:end);
    else 
        paddingFIFO(1:length(subbandSignal)) = subbandSignal;
        paddingFIFO = circshift(paddingFIFO,-length(subbandSignal));
    end
end