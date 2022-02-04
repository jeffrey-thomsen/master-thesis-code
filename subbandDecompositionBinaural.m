function [subbandSignals, AlgorithmParameters] = ...
  subbandDecompositionBinaural(inputSignal, AlgorithmParameters)

    analyzer.L = AlgorithmParameters.L.FilterStates.Gammatone.analyzer;
    analyzer.R = AlgorithmParameters.R.FilterStates.Gammatone.analyzer;
    
    [subbandSignals.L, analyzer.L] = ...
        subbandDecomposition(inputSignal(:,1), analyzer.L);

    [subbandSignals.R, analyzer.R] = ...
        subbandDecomposition(inputSignal(:,2), analyzer.R);

        
    AlgorithmParameters.L.FilterStates.Gammatone.analyzer = analyzer.L;
    AlgorithmParameters.R.FilterStates.Gammatone.analyzer = analyzer.R;

end