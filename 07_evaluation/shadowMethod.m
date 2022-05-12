% Shadow method: enhance target and interferer signals separately with 
% the SimulationData from the mixed signal processing for assessing
% SNR improvement
function [enhancedTarget, enhancedInterf] = shadowMethod(targetSignal, ...
  interfSignal, targetAngle, SimulationData, AlgorithmParameters, ...
  AlgorithmStates)

    AlgorithmParameters.targetRangeDeg = targetAngle + [-5 5];

    enhancedTarget = enhanceComparisonSignal('target', targetSignal, ...
        SimulationData, AlgorithmStates, AlgorithmParameters);

    enhancedInterf = enhanceComparisonSignal('interf', interfSignal, ...
        SimulationData, AlgorithmStates, AlgorithmParameters);

end