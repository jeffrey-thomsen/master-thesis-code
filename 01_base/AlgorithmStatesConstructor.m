function [AlgorithmStates, nBands] = AlgorithmStatesConstructor(AlgorithmParameters)
    
    States = struct;

    %% construct gammatone filterbanks
    [analyzer, synthesizer] = ...
        constructGammatoneFilterbank(AlgorithmParameters.Gammatone);
    
    GammatoneStates.analyzer = analyzer;
    GammatoneStates.synthesizer = synthesizer;

    nBands = length(GammatoneStates.analyzer.filters);

    %% preallocate states for all filters
    % seperate processing of each channel
    States.sigmaNormLP = 0;
    States.deltaNormLP = 0;
    States.ildLP = [0; 0];
    % combined binaural processing
    BinauralStates.ipdLP = 0;
    BinauralStates.ivsNumeratorLP = 0;
    BinauralStates.ivsDenominatorLP = 0;

    %% preallocate additional values that need to be updated continously

    nMaxSamplesP0Detection = ceil(AlgorithmParameters.Gammatone.samplingRateHz/...
        AlgorithmParameters.p0SearchRangeHz(1));
    % seperate processing of each channel
    States.p0DetectionFIFO = zeros(1,nMaxSamplesP0Detection);
    States.p0CandidateFIFO = zeros(nMaxSamplesP0Detection,1);
    States.instFreqPreviousValue = 0;
    % combined binaural processing
    BinauralStates.ivsPreviousValue = 0;
    %% generate separate set for each subband and channel

    L.GammatoneStates = GammatoneStates;
    R.GammatoneStates = GammatoneStates;

    [L.ProcessingStates{1:nBands}] = deal(States);
    [R.ProcessingStates{1:nBands}] = deal(States);
    [Binaural.ProcessingStates{1:nBands}] = deal(BinauralStates);

    AlgorithmStates.L = L;
    AlgorithmStates.R = R;
    AlgorithmStates.Binaural = Binaural;
end