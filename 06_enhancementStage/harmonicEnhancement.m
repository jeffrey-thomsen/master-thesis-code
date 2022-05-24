% Apply harmonic enhancement or cancellation to detected periodic subband
% samples
%
% Input:
% subbandSignal - array of subband signal samples
% iBand - index of subband to be processed
% ivsMask - cells containing logical arrays indicating the indices of the
% DOA estimates within each subband signal
% p0DetectedIndexVectors - cells containing information about position of 
% periodic samples within each subband signal and their respective detected
% period
% azimuthDegCells -  cells containing DOA estimate values for each subband
% snrDesired, sigmaDesired, deltaDesired - cells containing SNR, Sigma and
% Delta and SNR values of detected periodic subband samples
% AlgorithmParameters - struct containing parametric information for the
% simulation
%
% Output:
% subbandSignal - array of processed subband signal samples
% targetSampleIndices - vector containing indices of the subband samples
% that were identified as target speech and processed with harmonic
% enhancement
% interfererSampleIndices - vector containing indices of the subband 
% samples that were identified as interferer speech and processed with
% harmonic enhancement
function [subbandSignal, targetSampleIndices, interfererSampleIndices] = ...
  harmonicEnhancement(subbandSignal, iBand, ivsMask, ...
  p0DetectedIndexVectors, azimuthDegCells, snrDesired, sigmaDesired, ...
  deltaDesired, AlgorithmParameters)

    %% remove all samples (respective azimuth estimates) for which no 
    %% periodicity was detected

    if AlgorithmParameters.coherenceMask
        % create logical array true for samples that passed the IVS mask in
        % DOA estimation AND had periodicity detected
        coherentPeriodicComponentLogicalVector = ...
            ivsMask{iBand} & p0DetectedIndexVectors{iBand};
        % initialize array of length of ivsMask(=length of signal)
        azimuthDegVector = zeros(size(ivsMask{iBand}));
        % write azimuth values that passed IVS mask into array at correct
        % indices
        azimuthDegVector(ivsMask{iBand}) = azimuthDegCells{iBand};
        % extract azimuth values that satisfy both IVS and periodic condition
        coherentPeriodicAzimuthDegVector = ...
            azimuthDegVector(coherentPeriodicComponentLogicalVector);
    else
        % create logical array true for samples that had periodicity detected
        coherentPeriodicComponentLogicalVector = ...
            logical(p0DetectedIndexVectors{iBand});
        % extract azimuth values that satisfy periodic condition
        coherentPeriodicAzimuthDegVector = ...
            azimuthDegCells{iBand}(coherentPeriodicComponentLogicalVector);
    end
    

    %% evaluate target and interferer angle conditions
    targetConditionLogicalVectorReCoherentPeriodicSamples = ...
        coherentPeriodicAzimuthDegVector >= AlgorithmParameters.targetRangeDeg(1) & ...
        coherentPeriodicAzimuthDegVector <= AlgorithmParameters.targetRangeDeg(2);
    interfererConditionLogicalVectorReCoherentPeriodicSamples = ...
        coherentPeriodicAzimuthDegVector < AlgorithmParameters.targetRangeDeg(1) | ...
        coherentPeriodicAzimuthDegVector > AlgorithmParameters.targetRangeDeg(2);

    % extract signal index values for which IVS and periodic condition hold
    coherentPeriodicSampleIndices = ...
        find(coherentPeriodicComponentLogicalVector);
    % from that, extract signal index values for which target and 
    % interferer conditions hold
    targetSampleIndices = ...
        coherentPeriodicSampleIndices(targetConditionLogicalVectorReCoherentPeriodicSamples);
    interfererSampleIndices = ...
        coherentPeriodicSampleIndices(interfererConditionLogicalVectorReCoherentPeriodicSamples);
        
    % extract signal index values for which period was detected, i.e. sigma
    % and delta were computed
    periodicSampleIndices = find(p0DetectedIndexVectors{iBand});
    % compare those with target and interferer signal index samples to
    % create masks for extracting the correct sigma and delta values
    targetConditionLogicalVectorRePeriodicSamples = ...
        ismember(periodicSampleIndices, targetSampleIndices);
    interfererConditionLogicalVectorRePeriodicSamples = ...
        ismember(periodicSampleIndices, interfererSampleIndices);

    if AlgorithmParameters.bruemannSnrCondition
        % apply SNR condition according to Bruemann, applying a periodicity
        % threshold for cancellation but not enhancement
        assert(~AlgorithmParameters.snrCondition, 'Double SNR conditions (global and Bruemann) not allowed. Choose one!')
        snrConditionLogicalVectorRePeriodicSamples = ...
            snrDesired{iBand}>1;
        snrConditionPeriodicSampleIndices = ...
            periodicSampleIndices(snrConditionLogicalVectorRePeriodicSamples(1:numel(periodicSampleIndices)));
        
        interfererConditionLogicalVectorRePeriodicSamples = ...
            interfererConditionLogicalVectorRePeriodicSamples ...
            & snrConditionLogicalVectorRePeriodicSamples(1:numel(periodicSampleIndices))';
        interfererSampleIndices = ...
            interfererSampleIndices(ismember(interfererSampleIndices, snrConditionPeriodicSampleIndices));
    end

    
    if ~AlgorithmParameters.DOAProcessing
        % control condition: apply cancellation to all detected periodic
        % samples, regardless of coherence or angle of incidence
        interfererSampleIndices = periodicSampleIndices;
        interfererConditionLogicalVectorRePeriodicSamples = ...
            ismember(periodicSampleIndices, interfererSampleIndices);
        targetSampleIndices = [];
        targetConditionLogicalVectorRePeriodicSamples = [];
        
        % control condition: apply enhancement
        targetSampleIndices = periodicSampleIndices;
        targetConditionLogicalVectorRePeriodicSamples = ...
            ismember(periodicSampleIndices, targetSampleIndices);
        interfererSampleIndices = [];
        interfererConditionLogicalVectorRePeriodicSamples = [];
    end

    %% replace periodic samples from target and interferer angles with
    %% sigmas and deltas respectively

    if AlgorithmParameters.Enhancement
        subbandSignal(targetSampleIndices,iBand) = ...
            sigmaDesired{iBand}(targetConditionLogicalVectorRePeriodicSamples);
    end

    if AlgorithmParameters.Cancellation
        subbandSignal(interfererSampleIndices,iBand) = ...
            deltaDesired{iBand}(interfererConditionLogicalVectorRePeriodicSamples);
    end
end