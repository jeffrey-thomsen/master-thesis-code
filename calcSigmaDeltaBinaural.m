function [sigma, delta] = ...
  calcSigmaDeltaBinaural(subbandSignal, p0SearchRangeSamplesVector, ...
  mode, varargin)
    
    switch mode
        case 'range'
            [sigma.L, delta.L] = calcSigmaDeltaRange(subbandSignal.L, ...
                p0SearchRangeSamplesVector);
            [sigma.R, delta.R] = calcSigmaDeltaRange(subbandSignal.R, ...
                p0SearchRangeSamplesVector);
        case 'discrete'
            p0CandidateIndexVector = varargin{1};
            [sigma.L, delta.L] = calcSigmaDeltaDiscrete(subbandSignal.L,...
                p0SearchRangeSamplesVector, p0CandidateIndexVector.L);
            [sigma.R, delta.R] = calcSigmaDeltaDiscrete(subbandSignal.R,...
                p0SearchRangeSamplesVector, p0CandidateIndexVector.R);
        otherwise
            warning('must specify mode (range or discrete)')
    end


end

