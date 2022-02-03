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