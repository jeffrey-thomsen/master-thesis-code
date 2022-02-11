function [AlgorithmStates, nBands] = AlgorithmStatesConstructor(AlgorithmParameters)
    
    States = struct;

    %% construct gammatone filterbanks
    [analyzer, synthesizer] = ...
        constructGammatoneFilterbank(AlgorithmParameters.Gammatone);
    
    GammatoneStates.analyzer = analyzer;
    GammatoneStates.synthesizer = synthesizer;

    nBands = length(GammatoneStates.analyzer.filters);

    %% preallocate states for all filters
    States.sigmaNormLP = 0;
    States.deltaNormLP = 0;
    States.itfLP = 0;
    States.ildLP = 0;
    States.ivsNumeratorLP = 0;
    States.ivsDenominatorLP = 0;

    %% preallocate additional values that need to be updated continously

    nMaxSamplesP0Detection = ceil(AlgorithmParameters.Gammatone.samplingRateHz/...
        AlgorithmParameters.p0SearchRangeHz(1));

    States.p0DetectionFIFO = zeros(1,nMaxSamplesP0Detection);
    States.p0CandidateFIFO = zeros(1,nMaxSamplesP0Detection);
    States.ivsGradientPrevVal = 0;
    
    %% generate separate set for each subband and channel

    L.GammatoneStates = GammatoneStates;
    R.GammatoneStates = GammatoneStates;

    [L.ProcessingStates{1:nBands}] = deal(States);
    [R.ProcessingStates{1:nBands}] = deal(States);

    AlgorithmStates.L = L;
    AlgorithmStates.R = R;
end