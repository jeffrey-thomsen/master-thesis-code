function [sigma, delta] = ...
  calcSigmaDeltaRange(subbandSignal, p0SearchRangeSamplesVector)
    
    sigmaFunc=@(a,b) bsxfun(@plus, a,b.');
    deltaFunc=@(a,b) bsxfun(@minus,a,b.');
    pastFunc=@(s,p0Values) ...
        s(bsxfun(@minus,((1+max(p0Values)):end),p0Values)); %FIFO
    padding = zeros(max(p0SearchRangeSamplesVector) - ...
        min(p0SearchRangeSamplesVector)+1, max(p0SearchRangeSamplesVector));

    pastSamplesInSearchRange = ...
        [padding, pastFunc(subbandSignal, p0SearchRangeSamplesVector)];
    sigma = sigmaFunc(subbandSignal, pastSamplesInSearchRange);
    delta = deltaFunc(subbandSignal, pastSamplesInSearchRange);
end