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

    AlgorithmParameters.Gammatone = GammatoneParameters;

    AlgorithmParameters.p0SearchRangeHz = [100 500]; % Bruemann2018
    AlgorithmParameters.snrLPFilterTau = 0.04; % Bruemann2018
    AlgorithmParameters.ivsThreshold = 0.98; % Dietz2011
    AlgorithmParameters.nCyclesTau = 5; % Dietz2011
    AlgorithmParameters.snrThresholdInDb = 0;

    AlgorithmParameters.Cancellation = true;
    AlgorithmParameters.Enhancement = false;

end