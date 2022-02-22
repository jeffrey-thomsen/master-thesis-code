function p0CandidateIndexVector = subbandSnrPeakDetectionBinaural(snr)
    p0CandidateIndexVector.L = subbandSnrPeakDetection(snr.L);
    p0CandidateIndexVector.R = subbandSnrPeakDetection(snr.R);
end

function p0CandidateGlobalMaxIndexVector = subbandSnrPeakDetection(snr)
    % within snr matrix, find global maxima along p0 search range for each 
    % signal sample, excluding edge maxima
    % snr - NxM matrix - N: length of subband signal, M: number of values
    % in p0 search range
    % Returns vector of length N, containing for each signal sample: the 
    % index of maximum SNR value within p0 search range, or 0 if none was 
    % found for that sample


    % Calculate difference vector to find positive or negative gradient
    % assign +1 if positive gradient (difference)
    % assign -1 if negative gradient (difference)
    snrDifferenceGradient = sign(diff(snr,1,2)); % Set 1 for positive, -1 for negative gradient
    snrDifferenceGradient(snrDifferenceGradient==0) = 1; % if sign is 0, set to positive
    snrDiffereceGradientShifted = ...
        [snrDifferenceGradient(:,1), snrDifferenceGradient];
    snrDifferenceGradient = ...
        [snrDifferenceGradient, snrDifferenceGradient(:,end)];

    % Set local maxima to logical true if non-edge maximum detected
    p0CandidateLocalMaxLogicalMatrixIndices = false(size(snr));
    % maximum == sign-change from positive to negative
    p0CandidateLocalMaxLogicalMatrixIndices((snrDiffereceGradientShifted==1) & ...
        (snrDifferenceGradient==-1)) = true;

    % Pass on the amplitudes of the maximal SNRs at the relevant points, 
    % The rest should be set to 0.
    snrP0CandidateLocalMaxValues = zeros(size(snr));
    snrP0CandidateLocalMaxValues(p0CandidateLocalMaxLogicalMatrixIndices) = ...
        snr(p0CandidateLocalMaxLogicalMatrixIndices);

    % find maximum of SNR non-edge stationary points along FIFO arrays
    [~,p0CandidateGlobalMaxIndexVector] = ...
        max(snrP0CandidateLocalMaxValues,[],2);
    % If no maximum is found, the max() function returns the position 1 as
    % default, so this must be ensured to be removed again:
    p0CandidateGlobalMaxIndexVector(p0CandidateGlobalMaxIndexVector==1) = 0;
end