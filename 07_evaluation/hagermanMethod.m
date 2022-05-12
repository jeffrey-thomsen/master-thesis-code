% Hagerman ond Olofsson's (2004) method for calculating SNR improvement
function [enhancedTarget, enhancedInterf] = hagermanMethod(targetSignal, ...
  interfSignal, targetAngle, AlgorithmParameters, AlgorithmStates)
    
    testSignalPlus  = targetSignal + interfSignal;
    testSignalMinus = targetSignal - interfSignal;
    
    AlgorithmParameters.targetRangeDeg = targetAngle + [-5 5];
    
    [enhancedSignalPlus, ~, ~] = speechEnhancement(testSignalPlus, ...
        AlgorithmParameters, AlgorithmStates);
    
    [enhancedSignalMinus, ~, ~] = speechEnhancement(testSignalMinus, ...
        AlgorithmParameters, AlgorithmStates);
    
    enhancedTarget = (enhancedSignalPlus + enhancedSignalMinus)./2;
    enhancedInterf = (enhancedSignalPlus - enhancedSignalMinus)./2;
    
end
