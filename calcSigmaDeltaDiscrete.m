function [sigma, delta] = calcSigmaDeltaDiscrete(subbandSignal, ...
  p0SearchRangeSamplesVector, p0CandidateIndexVector)
    
    sigmaFunc=@(a,b) bsxfun(@plus, a,b);
    deltaFunc=@(a,b) bsxfun(@minus,a,b);
    pastFunc=@(s,p0Indices,p0Values) ...
        s(bsxfun(@minus,p0Indices,p0Values));

% at each p0DetectedIndex calculate sigma+delta with the p0DetectedSample
    p0CandidateSampleIndexVector = find(p0CandidateIndexVector); % indices of signal samples at which a p0 candidate was detected
    p0CandidateNonZeroIndexVector = p0CandidateIndexVector(p0CandidateIndexVector~=0); % indices of p0 candidate values 
    p0CandidateValueVector = p0SearchRangeSamplesVector(p0CandidateNonZeroIndexVector); % map indices to p0 candidate values (in samples)
    pastSamplesAtP0Candidate = pastFunc(subbandSignal, ...
        p0CandidateSampleIndexVector, p0CandidateValueVector);

    sigma = sigmaFunc(subbandSignal(p0CandidateSampleIndexVector), ...
        pastSamplesAtP0Candidate);
    delta = deltaFunc(subbandSignal(p0CandidateSampleIndexVector), ...
        pastSamplesAtP0Candidate);
end