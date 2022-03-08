% Apply enhancement/cancellation to identified samples within a binaural
% set of subband signals
%
% Input:
% subbandSignalArray - struct containing subband signal arrays of left and
% right channel
% azimuthDegCells - cells containing DOA estimate values for each subband
% ivsMask - cells containing logical arrays indicating the indices of the
% DOA estimates within each subband signal
% sigmaDesired, deltaDesired, snrDesired - structs containing Sigma, Delta
% and SNR values of detected periodic subband samples of left and right 
% channel
% p0DetectedIndexVectors - structs containing information about position of
% periodic samples within each subband signal and their respective detected
% period of left and right channel
% AlgorithmParameters - struct containing parametric information for the
% simulation
% 
% Output:
% enhancedSubbandSignalArray - struct containing the processed subband
% signal arrays of left and right channel
% targetSampleIndices - struct containing indices of the subband samples
% that were identified as target speech and processed with harmonic
% enhancement of left and right channel
% interfererSampleIndices - struct containing indices of the subband 
% samples that were identified as interferer speech and processed with
% harmonic cancellation of left and right channel
function [enhancedSubbandSignalArray, targetSampleIndices, ...
  interfererSampleIndices] = harmonicEnhancementBinaural(...
  subbandSignalArray, azimuthDegCells, ivsMask, sigmaDesired, deltaDesired,...
  snrDesired, p0DetectedIndexVectors, AlgorithmParameters)

    nBands = AlgorithmParameters.Gammatone.nBands;
    for iBand = 1:nBands

        [subbandSignalArray.L, targetSampleIndices.L{iBand}, ...
            interfererSampleIndices.L{iBand}] = harmonicEnhancement(...
            subbandSignalArray.L, iBand, ivsMask, ...
            p0DetectedIndexVectors.L, azimuthDegCells, snrDesired.L, ...
            sigmaDesired.L, deltaDesired.L, AlgorithmParameters);
        
        [subbandSignalArray.R, targetSampleIndices.R{iBand}, ...
            interfererSampleIndices.R{iBand}] = harmonicEnhancement(...
            subbandSignalArray.R, iBand, ivsMask, ...
            p0DetectedIndexVectors.R, azimuthDegCells, snrDesired.R, ...
            sigmaDesired.R, deltaDesired.R, AlgorithmParameters);

% conditions from Bruemann2018
%% If maxSNR below threshold and angle within range, pass on sigmas 
%% If maxSNR above threshold and angle within range, pass on sigmas
%% If maxSNR above threshold and angle outside of range, pass on deltas

    end

    enhancedSubbandSignalArray = subbandSignalArray;
end