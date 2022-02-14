function [sigma, delta, FIFO] = calcSigmaDeltaRange(subbandSignal, ...
  p0SearchRangeSamplesVector, FIFO)
    % Computes sigmas and deltas of subbandSignal for a range of previous 
    % signal samples stored in FIFO. 
    % FIFO updated every sample.
    % Range of past samples for sigma and delta computation specified in
    % p0SearchRangeSamplesVector.
    % Returns arrays sigma and delta of shape NxM - N: length of signal, M:
    % number of samples in p0 search range
    nSigSamples = length(subbandSignal);
    nP0Values = length(p0SearchRangeSamplesVector);
    sigma = zeros(nSigSamples,nP0Values);
    delta = zeros(nSigSamples,nP0Values);
    for i = 1:nSigSamples
        sigma(i,:) = subbandSignal(i) + ...
            FIFO(p0SearchRangeSamplesVector);
        delta(i,:) = subbandSignal(i) - ...
            FIFO(p0SearchRangeSamplesVector);
        % update FIFO
        FIFO(end) = subbandSignal(i);
        FIFO = circshift(FIFO,1);
    end
end


