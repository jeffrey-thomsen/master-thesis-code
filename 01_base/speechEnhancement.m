% run a binaural signal, signal block or signal sample through the
% Thomsen2022 speech enhancement algorithm
%
% Inpput:
% inputSignal - real-valued Nx2 matrix which must satisfy the conditions 
% that can be checked by running it through testInputSignal.m
% AlgorithmParameters - struct containing all simulation parameters
% AlgorithmStates - struct containing filter states and FIFO arrays that
% need to be handed over to the next signal block/sample to be processed
%
% Output:
% enhancedSignal - real-valued Nx2 matrix containing the processed binaural
% signal
% AlgorithmStates - see above
% simulationData - struct containing intermediate results from the
% different algorithm stages for evaluating its performance
function [enhancedSignal, AlgorithmStates, SimulationData] = ...
  speechEnhancement(inputSignal, AlgorithmParameters, AlgorithmStates)

    % Gammatone analysis filterbank - decompose signal into frequency bands
    [subbandSignalArray, AlgorithmStates] = ...
        subbandDecompositionBinaural(inputSignal, AlgorithmStates);
    
    if AlgorithmParameters.Cancellation || AlgorithmParameters.Enhancement

        % Periodicity analysis - detect samples with periodic components
        [sigmaDesired, deltaDesired, snrDesired, p0DetectedIndexVectors, ...
            AlgorithmStates, p0SearchRangeSamplesVector] = ...
            periodicityAnalysis(subbandSignalArray, ...
            AlgorithmParameters, AlgorithmStates);
        
        % DOA estimation - estimate angle of incidence for coherent sources
        [ipdRadCells, ivsMaskCells, ~, AlgorithmStates, ...
            ipdDisambiguatedLogicalCells, itdSecCells, ...
            itdDisambiguatedLogicalCells] = ...
            interauralFeatureComputation(subbandSignalArray, ...
            AlgorithmParameters, AlgorithmStates);

        azimuthDegCells = azimuthEstimation(ipdRadCells, itdSecCells, ...
            ivsMaskCells, AlgorithmParameters);
        
        % signal enhancement - 
        % enhance periodic samples in target angle and
        % suppress periodic samples outside of target angle
        [enhancedSubbandSignalArray, targetSampleIndices, ...
            interfSampleIndices] = harmonicEnhancementBinaural(...
            subbandSignalArray, azimuthDegCells, ivsMaskCells, ...
            sigmaDesired, deltaDesired, snrDesired, ...
            p0DetectedIndexVectors, AlgorithmParameters);
        
        % store data for evaluating simulation
        SimulationData.p0DetectedIndexVectors = p0DetectedIndexVectors;
        SimulationData.p0SearchRangeSamplesVector = p0SearchRangeSamplesVector;
        SimulationData.ipdRadCells = ipdRadCells;
        SimulationData.ivsMaskCells = ivsMaskCells;
        SimulationData.ipdDisambiguatedLogicalCells = ipdDisambiguatedLogicalCells;
        SimulationData.azimuthDegCells = azimuthDegCells;
        SimulationData.targetSampleIndices = targetSampleIndices;
        SimulationData.interfSampleIndices = interfSampleIndices;
        SimulationData.itdSecCells = itdSecCells;
        SimulationData.itdDisambiguatedLogicalCells = itdDisambiguatedLogicalCells;
    
        SimulationData.sigmaDesired = sigmaDesired;
        SimulationData.deltaDesired = deltaDesired;
        SimulationData.snrDesired = snrDesired;

    else % just for testing the gammatone filterbank without processing
        enhancedSubbandSignalArray = subbandSignalArray;
        SimulationData.msg = 'No simulation data. Only gammatone decomposition and resynthesis.';
    end
    
    % Gammatone synthesis filterbank - resynthesize enhanced subband signals
    [enhancedSignal, AlgorithmStates] = ...
        subbandResynthesisBinaural(enhancedSubbandSignalArray, ...
        AlgorithmStates);

   
end