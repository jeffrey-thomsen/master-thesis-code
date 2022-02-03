function [subbandSignals, analyzer] = ...
  subbandDecomposition(inputSignal, analyzer)

    [subbandSignals, analyzer] = ...
        hohmann2002_process(analyzer, inputSignal);

end