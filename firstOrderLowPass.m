% This is a simple low-pass filter y(n) = (1-a)*x(n) - a*y(n-1)
% Meaning of parameter a:
% a - damping coefficient (0 - no filtering, 1 - flat output)
%
% tau - time constant in seconds, when the filter decays to exp(-1)
% a = exp(-1/(fs*tau)) where fs - sampling frequency, or identical as used
% here: a = exp(-(T/tau)) where T = 1/fs - sampling period in seconds  
%
function [outputSignal, filterState] = firstOrderLowPass(inputSignal, ...
  filterState, AlgorithmParameters, varargin)

    %% Filter variable definition
    if ~isempty(varargin)
        tau = varargin{1};
    else
        tau = AlgorithmParameters.snrLPFilterTau;
    end
    fs = AlgorithmParameters.Gammatone.samplingRateHz;
    T = 1/fs;
    a = exp(-(T/tau)); % filter constant a: (unit-less (seconds/seconds))
    %% Preallocation of output signal
    outputSignal = zeros(size(inputSignal)); 
    %% Filter iterations
    outputSignal(1,:) = a*filterState + (1-a)*inputSignal(1,:);
    % rest of averaged elements (using past averaged elements)
    for iSample = 2:size(outputSignal,1)
        outputSignal(iSample,:) = ...
            a*outputSignal(iSample-1,:) + (1-a)*inputSignal(iSample,:);
    end
    filterState = outputSignal(end,:);
end