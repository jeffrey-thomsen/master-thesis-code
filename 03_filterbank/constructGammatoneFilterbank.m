% Generates gammatone analysis and synthesis filterbank according to
% Hohmann2002
% Parameters - struct containing all parameters needed to define the
% filterbank
% analyzer - struct containing parameters and states of the analysis
% filterbankd
% synthesizer - struct containing parameters and states of the synthesis
% filterbankd
function [analyzer, synthesizer] = constructGammatoneFilterbank(Parameters)

    fLow = Parameters.fLow;
    fHigh = Parameters.fHigh;
    baseFreqHz = Parameters.baseFreqHz;
    samplingRateHz = Parameters.samplingRateHz;
    filtersPerErb = Parameters.filtersPerErb;
    desiredDelayInSeconds = Parameters.desiredDelayInSeconds;
    filterOrder = Parameters.filterOrder;
    bandwidthFactor = Parameters.bandwidthFactor;

    analyzer = hohmann2002(samplingRateHz, fLow, baseFreqHz, ...
        fHigh, filtersPerErb, filterOrder, bandwidthFactor);

    synthesizer = hohmann2002_synth(analyzer, desiredDelayInSeconds);
    
end