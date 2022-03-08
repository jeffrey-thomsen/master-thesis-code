% Calculate sigma and delta values from a set of binaural subband signals
% for period values specified in p0SearchRangeSamplesVector
%
% Input:
% subbandSignal - struct containing a subband signal from the left and
% right channel
% p0SearchRangeSamplesVector - index range vector specifying the p0 search
% range period values in samples
% FilterStates - struct containing the FIFOs used to store the past signal
% samples needed for sigma and delta calculation
% mode - string 'range' or 'discrete' specifying the operation mode
% varargin - if mode is 'discrete', must be a struct containing p0
% candidate index vectors for the left and right channel, specifiyng which
% period to apply to which subband samples in discrete sigma and delta
% calculation
%
% Output:
% FilterStates - see above
% delta, sigma - structs containing the calculated values for the left and
% right channel
function [sigma, delta, FilterStates] = ...
  calcSigmaDeltaBinaural(subbandSignal, p0SearchRangeSamplesVector, ...
  FilterStates, mode, varargin)
    
    switch mode
        case 'range'
            [sigma.L, delta.L, FilterStates.L.p0DetectionFIFO] = ...
                calcSigmaDeltaRange(subbandSignal.L, ...
                p0SearchRangeSamplesVector, FilterStates.L.p0DetectionFIFO);
            [sigma.R, delta.R, FilterStates.R.p0DetectionFIFO] = ...
                calcSigmaDeltaRange(subbandSignal.R, ...
                p0SearchRangeSamplesVector, FilterStates.R.p0DetectionFIFO);
        case 'discrete'
            p0CandidateIndexVector = varargin{1};
            [sigma.L, delta.L, FilterStates.L.p0CandidateFIFO] = ...
                calcSigmaDeltaDiscrete(subbandSignal.L,...
                p0SearchRangeSamplesVector, p0CandidateIndexVector.L,...
                FilterStates.L.p0CandidateFIFO);
            [sigma.R, delta.R, FilterStates.R.p0CandidateFIFO] = ...
                calcSigmaDeltaDiscrete(subbandSignal.R,...
                p0SearchRangeSamplesVector, p0CandidateIndexVector.R,...
                FilterStates.R.p0CandidateFIFO);
        otherwise
            warning('must specify mode (range or discrete)')
    end


end

