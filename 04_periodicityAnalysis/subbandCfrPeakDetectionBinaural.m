% Find the maximum CFR (periodicity) value among a range of period values
% for each signal sample
%
% Input:
% cfr - struct containing arrays of CFR values of the left and right
% channel
%
% Output:
% p0CandidateVector - struct containing the indices of the maximal CFR
% value for each subband signal sample, or 0 if no maximum is found for a
% sample, of the left and right channel
function p0CandidateIndexVector = subbandCfrPeakDetectionBinaural(cfr, ...
  AlgorithmParameters)

    p0CandidateIndexVector.L = ...
        subbandCfrPeakDetection(cfr.L, AlgorithmParameters);
    p0CandidateIndexVector.R = ...
        subbandCfrPeakDetection(cfr.R, AlgorithmParameters);
end

function p0CandidateGlobalMaxIndexVector = ...
  subbandCfrPeakDetection(cfr, AlgorithmParameters)
    % within cfr matrix, find global maxima along p0 search range for each 
    % signal sample, excluding edge maxima
    % cfr - NxM real-valued matrix - N: length of subband signal, M: number
    % of values in p0 search range
    % Returns vector of length N, containing for each signal sample: the 
    % index of maximum CFR value within p0 search range, or 0 if none was 
    % found for that sample

    if ~isreal(cfr)
        warning('CFR values for subband CFR peak detection must be real-valued!')
    end
    % Calculate difference vector to find positive or negative gradient
    % assign +1 if positive gradient (difference)
    % assign -1 if negative gradient (difference)
    cfrDifferenceGradient = sign(diff(cfr,1,2)); % Set 1 for positive, -1 for negative gradient
    cfrDifferenceGradient(cfrDifferenceGradient==0) = 1; % if sign is 0, set to positive
    cfrDiffereceGradientShifted = ...
        [cfrDifferenceGradient(:,1), cfrDifferenceGradient];
    cfrDifferenceGradient = ...
        [cfrDifferenceGradient, cfrDifferenceGradient(:,end)];

    % Set local maxima to logical true if non-edge maximum detected
    p0CandidateLocalMaxLogicalMatrixIndices = false(size(cfr));
    % maximum == sign-change from positive to negative
    p0CandidateLocalMaxLogicalMatrixIndices((cfrDiffereceGradientShifted==1) & ...
        (cfrDifferenceGradient==-1)) = true;

    % Pass on the amplitudes of the maximal CFRs at the relevant points, 
    % The rest should be set to 0.
    cfrP0CandidateLocalMaxValues = zeros(size(cfr));
    cfrP0CandidateLocalMaxValues(p0CandidateLocalMaxLogicalMatrixIndices) = ...
        cfr(p0CandidateLocalMaxLogicalMatrixIndices);
    
    if AlgorithmParameters.cfr0mask
        % convert CFR0 threshold from dB to 1 scale
        cfr0Threshold = 10^(0.1*AlgorithmParameters.cfr0ThresholdInDb);
        % Apply CFR0 filter mask
        cfrP0CandidateLocalMaxValues(...
            cfrP0CandidateLocalMaxValues < cfr0Threshold) = 0;
    end

    % find maximum of CFR non-edge stationary points along FIFO arrays
    [~,p0CandidateGlobalMaxIndexVector] = ...
        max(cfrP0CandidateLocalMaxValues,[],2);
    % If no maximum is found, the max() function returns the position 1 as
    % default, so this must be ensured to be removed again:
    p0CandidateGlobalMaxIndexVector(p0CandidateGlobalMaxIndexVector==1) = 0;

    % control condition for evaluation
    if AlgorithmParameters.RandomP0
        p0CandidateGlobalMaxIndexVector(p0CandidateGlobalMaxIndexVector~=0) = randi(size(cfr,2),nnz(p0CandidateGlobalMaxIndexVector),1);
    end
end