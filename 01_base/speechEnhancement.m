% run a binaural signal, signal block or signal sample through the
% Thomsen2022 speech enhancement algorithm
% inputSignal - real-valued Nx2 matrix which must satisfy the conditions 
% that can be checked by running it through testInputSignal.m
% AlgorithmParameters - struct containing all information necessary for
% running the simulation including filter states
% enhancedSignal - real-valued Nx2 matrix containing the processed binaural
% signal
function [enhancedSignal, AlgorithmStates, simulationData] = ...
  speechEnhancement(inputSignal, AlgorithmParameters, AlgorithmStates)

    % Gammatone analysis filterbank - decompose signal into frequency bands
    [subbandSignalArray, AlgorithmStates] = ...
        subbandDecompositionBinaural(inputSignal, AlgorithmStates);
    
    if AlgorithmParameters.Cancellation || AlgorithmParameters.Enhancement
        % Periodicity analysis - detect subband samples with periodic components
        [sigmaDesired, deltaDesired, snrDesired, p0DetectedIndexVectors, ...
            AlgorithmStates, p0SearchRangeSamplesVector] = periodicityAnalysis(subbandSignalArray, ...
            AlgorithmParameters, AlgorithmStates);
        
        % DOA estimation - estimate angle of incidence for coherent sources
        [ipdRadCells, ivsMaskCells, ~, AlgorithmStates, ...
            ipdDisambiguatedLogicalCells, itdSecCells, ...
            itdDisambiguatedLogicalCells] = ...
            interauralFeatureComputation(subbandSignalArray, ...
            AlgorithmParameters, AlgorithmStates);

        azimuthDegCells = azimuthEstimation(ipdRadCells, itdSecCells, ...
            ivsMaskCells, AlgorithmParameters);
        
        % signal enhancement - apply harmonic cancellation to unwanted 
        % periodic components
        [enhancedSubbandSignalArray, targetSampleIndices, interfSampleIndices] = ...
            harmonicEnhancementBinaural(subbandSignalArray, azimuthDegCells, ...
            ivsMaskCells, sigmaDesired, deltaDesired, snrDesired, ...
            p0DetectedIndexVectors, AlgorithmParameters);
        
        % store data for evaluating simulation
        simulationData.p0DetectedIndexVectors = p0DetectedIndexVectors;
        simulationData.p0DetectedIndexVectors = p0DetectedIndexVectors;
        simulationData.p0SearchRangeSamplesVector = p0SearchRangeSamplesVector;
        simulationData.ipdRadCells = ipdRadCells;
        simulationData.ivsMaskCells = ivsMaskCells;
        simulationData.ipdDisambiguatedLogicalCells = ipdDisambiguatedLogicalCells;
        simulationData.azimuthDegCells = azimuthDegCells;
        simulationData.targetSampleIndices = targetSampleIndices;
        simulationData.interfSampleIndices = interfSampleIndices;
        simulationData.itdSecCells = itdSecCells;
        simulationData.itdDisambiguatedLogicalCells = itdDisambiguatedLogicalCells;
    
    else % just for testing the gammatone filterbank without processing
        enhancedSubbandSignalArray = subbandSignalArray;
    end
    
    % Gammatone synthesis filterbank - resynthesize enhanced subband signals
    [enhancedSignal, AlgorithmStates] = ...
        subbandResynthesisBinaural(enhancedSubbandSignalArray, ...
        AlgorithmStates);

   
end