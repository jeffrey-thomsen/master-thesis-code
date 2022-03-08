% Runs the Hohmann2002 gammatone analysis filterbank
%
% Input:
% inputSignal - Nx1 real-valued signal
% analyzer - struct containing the filterbank parameters and states
%
% Output:
% subbandSignals - NxM complex-valued set of frequency subband signals
function [subbandSignals, analyzer] = ...
  subbandDecomposition(inputSignal, analyzer)

    [subbandSignals, analyzer] = ...
        hohmann2002_process(analyzer, inputSignal);

end