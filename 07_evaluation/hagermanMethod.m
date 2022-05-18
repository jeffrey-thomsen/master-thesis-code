% Hagerman ond Olofsson's (2004) method for calculating SNR improvement
function [enhancedTarget, enhancedInterf] = hagermanMethod(processingMode, ...
    targetAngle, AlgorithmParameters, AlgorithmStates, varargin)
    
    AlgorithmParameters.targetRangeDeg = targetAngle + [-5 5];

    switch processingMode
        case 'full-calc' % construct test signals and process both

            testSignalPlus  = varargin{1} + varargin{2};
            testSignalMinus = varargin{1} - varargin{2};
            
            [enhancedSignalPlus, ~, ~] = speechEnhancement(testSignalPlus, ...
                AlgorithmParameters, AlgorithmStates);
            
            [enhancedSignalMinus, ~, ~] = speechEnhancement(testSignalMinus, ...
                AlgorithmParameters, AlgorithmStates);

        case 'pre-calc' % receive pre-calculated regular enhanced signal, 
            % pre-constructed Hagerman test signal, saving one redundant
            % processing step
            % Note: This has been informally validated as equal to the 
            % 'full-calc' method above within numerical precision
            enhancedSignalPlus = varargin{1};

            [enhancedSignalMinus, ~, ~] = speechEnhancement(varargin{2}, ...
                AlgorithmParameters, AlgorithmStates);
    end

    enhancedTarget = (enhancedSignalPlus + enhancedSignalMinus)./2;
    enhancedInterf = (enhancedSignalPlus - enhancedSignalMinus)./2;
    
end
