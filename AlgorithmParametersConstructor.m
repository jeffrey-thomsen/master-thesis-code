function AlgorithmParameters = AlgorithmParametersConstructor()
    AlgorithmParameters = struct;
    
    GammatoneParameters.fLow = 70;
    GammatoneParameters.fHigh = 6700;
    GammatoneParameters.baseFreqHz = 1000;
    GammatoneParameters.samplingRateHz = 16000;
    GammatoneParameters.filtersPerErb = 1;
    GammatoneParameters.desiredDelayInSeconds = 0.004; % 0.016 ensures perfect reconstruction according to Chen2015
    GammatoneParameters.filterOrder = 4;
    GammatoneParameters.bandwidthFactor = 1.0;

    AlgorithmParameters.GammatoneParameters = GammatoneParameters;

    AlgorithmParameters.p0SearchRangeHz = [100 500]; % Bruemann2018
    AlgorithmParameters.snrLPFilterTau = 0.04; % Bruemann2018
    AlgorithmParameters.ivsThreshold = 0.98; % Dietz2011
    AlgorithmParameters.snrThresholdInDb = 0;

    AlgorithmParameters.Cancellation = true;
    AlgorithmParameters.Enhancement = false;
    
    %% construct gammatone filterbanks
    [analyzer, synthesizer] = ...
        constructGammatoneFilterbank(GammatoneParameters);
    
    FilterStates.Gammatone.analyzer = analyzer;
    FilterStates.Gammatone.synthesizer = synthesizer;

    AlgorithmParameters.GammatoneParameters.nBands = ...
        length(FilterStates.Gammatone.analyzer.filters);

    %% preallocate states for all filters
    FilterStates.sigmaNormLP = 0;
    FilterStates.deltaNormLP = 0;
    FilterStates.itfLP = 0;
    FilterStates.ildLP = 0;
    FilterStates.ivsNumeratorLP = 0;
    FilterStates.ivsDenominatorLP = 0;
    
    %% generate separate set for each channel
    AlgorithmParameters.L.FilterStates = FilterStates;
    AlgorithmParameters.R.FilterStates = FilterStates;

    AlgorithmParameters.L.p0DetectionFIFO = 0;
    AlgorithmParameters.R.p0DetectionFIFO = 0;

    AlgorithmParameters.L.ivsGradientPrevVal = 0;
    AlgorithmParameters.R.ivsGradientPrevVal = 0;

    
end