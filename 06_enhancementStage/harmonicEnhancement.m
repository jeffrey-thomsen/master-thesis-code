function [subbandSignal, targetSampleIndices, interfererSampleIndices] = harmonicEnhancement(snrDesired,...
  ivsMask, p0DetectedIndexVectors, azimuthDegCells, subbandSignal,...
  sigmaDesired, deltaDesired, iBand, AlgorithmParameters)

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
    targetConditionLogicalVectorReCoherentPeriodicSamples = ... % target in front
        abs(coherentPeriodicAzimuthDegVector) <= 5;
    interfererConditionLogicalVectorReCoherentPeriodicSamples = ... % interference from the sides
        abs(coherentPeriodicAzimuthDegVector) > 5;
%         coherentPeriodicAzimuthDegVector > 50 & ...
%         coherentPeriodicAzimuthDegVector < 70;

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

    if AlgorithmParameters.snrCondition
        % apply additional SNR condition
        snrConditionLogicalVectorRePeriodicSamples = ...
            snrDesired{iBand}>1; % AlgorithmParameters.snrThresholdInDb
        snrConditionPeriodicSampleIndices = ...
            periodicSampleIndices(snrConditionLogicalVectorRePeriodicSamples);
        
        interfererConditionLogicalVectorRePeriodicSamples = ...
            interfererConditionLogicalVectorRePeriodicSamples ...
            & snrConditionLogicalVectorRePeriodicSamples;
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