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

    AlgorithmParameters.lookuptableType = 'itd'; % 'ipd' or 'itd' to azimuth mapping

    AlgorithmParameters.p0SearchRangeHz = [100 350]; % Bruemann2018 100 500
    AlgorithmParameters.cfrLPFilterTau = 0.04; % Bruemann2018
    AlgorithmParameters.ivsThreshold = 0.98; % Dietz2011
    AlgorithmParameters.nCyclesTau = 5; % Dietz2011
    AlgorithmParameters.cfr0ThresholdInDb = 10; % Bruemann2018 0 dB
    
    AlgorithmParameters.targetRangeDeg = [-5 5]; % initialized at 0°

    AlgorithmParameters.ChenP0Detection = false;
    AlgorithmParameters.coherenceMask = true;
    AlgorithmParameters.azimuthPooling = false;
    AlgorithmParameters.cfr0mask = true;
    AlgorithmParameters.DOAProcessing = true;
    AlgorithmParameters.bruemannSimpleDOA = false;
    AlgorithmParameters.Cancellation = true;
    AlgorithmParameters.Enhancement = true;

    AlgorithmParameters.RandomP0 = false;


end