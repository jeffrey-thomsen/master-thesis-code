function [enhancedSubbandSignalArray, targetSampleIndices, interfererSampleIndices] = ...
  harmonicEnhancementBinaural(subbandSignalArray, azimuthDegCells, ivsMask, ...
  sigmaDesired, deltaDesired, snrDesired, p0DetectedIndexVectors,...
  AlgorithmParameters)

    nBands = AlgorithmParameters.Gammatone.nBands;
    for iBand = 1:nBands

        [subbandSignalArray.L, targetSampleIndices.L{iBand}, interfererSampleIndices.L{iBand}] = harmonicEnhancement(snrDesired.L,...
            ivsMask, p0DetectedIndexVectors.L, azimuthDegCells, ...
            subbandSignalArray.L, sigmaDesired.L, deltaDesired.L, iBand,...
            AlgorithmParameters);
        
        [subbandSignalArray.R, targetSampleIndices.R{iBand}, interfererSampleIndices.R{iBand}] = harmonicEnhancement(snrDesired.R,...
            ivsMask, p0DetectedIndexVectors.R, azimuthDegCells, ...
            subbandSignalArray.R, sigmaDesired.R, deltaDesired.R, iBand,...
            AlgorithmParameters);

% If SNR below specified level: logical 0, else: 1

%% If maxSNR below threshold and angle within range, pass on sigmas 
%% If maxSNR above threshold and angle within range, pass on sigmas
%% If maxSNR above threshold and angle outside of range, pass on deltas

    end

    enhancedSubbandSignalArray = subbandSignalArray;
end