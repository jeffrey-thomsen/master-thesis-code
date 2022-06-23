% Computes sigmas and deltas for specific samples of the subbandSignal,
% each for a specific p0 period value.
% FIFO updated at the end of processing for handing it to the next signal
% block/sample.
% Specified in p0CandidateIndexVector are the indices pointing to that
% value in the past within the p0SearchRangeSamplesVector, stored at
% the position corresponding to the signal sample index of
% subbandSignal for which the specific sigma and delta are computed.
function [sigma, delta, paddingFIFO] = calcSigmaDeltaDiscrete(subbandSignal, ...
  p0SearchRangeSamplesVector, p0CandidateIndexVector, paddingFIFO)

    sigmaFunc=@(a,b) bsxfun(@plus, a,b);
    deltaFunc=@(a,b) bsxfun(@minus,a,b);
    pastFunc=@(signal,p0IndicesReSignal,p0Values,shift) ...
        signal(bsxfun(@minus,p0IndicesReSignal+shift,p0Values));
    
    paddingFIFOshift = length(paddingFIFO);

    % indices of signal samples at which periodicity was detected
    p0CandidateIndexVectorReSignalSamples = find(p0CandidateIndexVector);
    % indices of corresponding p0 candidate values 
    p0CandidateIndexVectorReP0SearchRange = ...
        nonzeros(p0CandidateIndexVector);

    % map indices to p0 candidate values (in samples)
    p0CandidateValueVector = p0SearchRangeSamplesVector(p0CandidateIndexVectorReP0SearchRange);
    
    % read signal samples corresponding to the detected period p0 in the
    % past
    paddedSubbandSignal = [paddingFIFO; subbandSignal];
    pastSamplesAtP0Candidate = pastFunc(paddedSubbandSignal, ...
        p0CandidateIndexVectorReSignalSamples, p0CandidateValueVector, paddingFIFOshift);

    % calculate discrete sigma and delta values from the past signal
    % samples
    sigma = sigmaFunc(subbandSignal(p0CandidateIndexVectorReSignalSamples), ...
        pastSamplesAtP0Candidate);
    delta = deltaFunc(subbandSignal(p0CandidateIndexVectorReSignalSamples), ...
        pastSamplesAtP0Candidate);

    % update padding FIFO
    if paddingFIFOshift<length(subbandSignal)
        paddingFIFO = subbandSignal(end-paddingFIFOshift+1:end);
    else 
        paddingFIFO(1:length(subbandSignal)) = subbandSignal;
        paddingFIFO = circshift(paddingFIFO,-length(subbandSignal));
    end
end