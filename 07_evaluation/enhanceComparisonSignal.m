% Use previously computed information about time-frequency
% target/interferer glimpses, stored in SimulationData, to process a
% different binaural signal (of equal length) for comparison
function enhancedSignal = enhanceComparisonSignal(type, inputSignal, ...
    SimulationData, AlgorithmStates, AlgorithmParameters)

switch type
    case 'target'
        AlgorithmParameters.Cancellation    = false;
        AlgorithmParameters.Enhancement     = true;
    case 'interf'
        AlgorithmParameters.Cancellation    = true;
        AlgorithmParameters.Enhancement     = false;
end

subbandSignal = subbandDecompositionBinaural(inputSignal, AlgorithmStates);

nBands = AlgorithmParameters.Gammatone.nBands;
for iBand = 1:nBands
    States.L = AlgorithmStates.L.ProcessingStates{iBand};
    States.R = AlgorithmStates.R.ProcessingStates{iBand};
    p0Candidates.L = SimulationData.p0DetectedIndexVectors.L{iBand};
    p0Candidates.R = SimulationData.p0DetectedIndexVectors.R{iBand};
    tmpSubbandSignal.L = subbandSignal.L(:,iBand);
    tmpSubbandSignal.R = subbandSignal.R(:,iBand);
    [p0CandidateSigma, p0CandidateDelta] = ...
        calcSigmaDeltaBinaural(tmpSubbandSignal, ...
        SimulationData.p0SearchRangeSamplesVector, States, 'discrete', ...
        p0Candidates);
    sigmaCells.L{iBand} = p0CandidateSigma.L;
    sigmaCells.R{iBand} = p0CandidateSigma.R;
    deltaCells.L{iBand} = p0CandidateDelta.L;
    deltaCells.R{iBand} = p0CandidateDelta.R;
end

enhancedSubbandSignalArray = harmonicEnhancementBinaural(...
            subbandSignal, ...
            SimulationData.azimuthDegCells, ...
            SimulationData.ivsMaskCells, ...
            sigmaCells, ...
            deltaCells, ...
            SimulationData.cfrDesired, ...
            SimulationData.p0DetectedIndexVectors, ...
            AlgorithmParameters);

enhancedSignal = subbandResynthesisBinaural(enhancedSubbandSignalArray, ...
    AlgorithmStates);

end