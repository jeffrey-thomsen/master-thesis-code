function [outputSignal, filterState] = firstOrderLowPass(inputSignal, ...
  filterState, AlgorithmParameters)

    %% Filter variable definition
    tau = AlgorithmParameters.snrLPFilterTau;
    fs = AlgorithmParameters.Gammatone.samplingRateHz;
    T = 1/fs;
    a = exp(-(T/tau)); % filter constant a: (unit-less (seconds/seconds))
    %% Preallocation of output signal
    outputSignal = zeros(size(inputSignal)); 
    %% Filter iterations
    outputSignal(1,:) = a*filterState + (1-a)*inputSignal(1,:);
    % rest of averaged elements (using past averaged elements)
    for iSample = 2:length(outputSignal)
        outputSignal(iSample,:) = ...
            a*outputSignal(iSample-1,:) + (1-a)*inputSignal(iSample,:);
    end
    filterState = outputSignal(end,:);
end