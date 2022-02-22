% generate an IPD-to-azimuth mapping function in form of a polynomial
% hrtf - SOFA file
% AlgorithmParameters - struct
function [lookuptable, ivsArray] = ...
  ipdToAzimuthLookuptable(hrtf, AlgorithmParameters)

    %% Configuration
    
    samplingRateHz = AlgorithmParameters.Gammatone.samplingRateHz;
    
    % length of noise signal used for the calculation (samples)
    nSamples = 1*samplingRateHz;
    
    % initialize filter states for interaural feature computation
    [FilterStates, nBands] = ...
        AlgorithmStatesConstructor(AlgorithmParameters);
    AlgorithmParameters.Gammatone.nBands = nBands;
    
    Obj = hrtf;
    
    polynomialOrder = 12; % Dietz2011 paper says 9;
    
    %% Preparation
    
    % generate noise signal
    noiseSignal = randn(nSamples,1);
    % normalize noise signal
    noiseSignal = noiseSignal - mean(noiseSignal);
    noiseSignal = noiseSignal/max(abs(noiseSignal));
    % alternative tonal signal
    % noiseSignal = sin(2*pi*10000*(1/samplingRateHz:1/samplingRateHz:1))';
    
    % use only the -90 to 90 degree part of the HRIR set
    elevationAngleDeg = 0;
    idFrontHor = ...
        Obj.SourcePosition(:,2) == elevationAngleDeg & ... % horizontal plane
        (((Obj.SourcePosition(:,1) >= -90 & Obj.SourcePosition(:,1) < 0) ... % front right
        | Obj.SourcePosition(:,1) >= 270) | ... % alternative front right
        (Obj.SourcePosition(:,1) >= 0 & Obj.SourcePosition(:,1) <= 90)); % front left
    azimuthAngleDeg = Obj.SourcePosition(idFrontHor,1);
    % translate angles into -90...90 degree range
    azimuthAngleDeg(azimuthAngleDeg>180) = ...
        azimuthAngleDeg(azimuthAngleDeg>180)-360;
    
    %% Calculation 
    
    nAzimuth = length(azimuthAngleDeg);
    
    ipdRadMedianArray = cell(1,nBands); %zeros(nAzimuth, nBands);
%     ipdIvsMedianArray = cell(1,nBands); %zeros(nAzimuth, nBands);
    ivsCellsCells = cell(1, nAzimuth);
    
    for iAzimuth = 1:nAzimuth
    
        % generate noise coming from the given direction
        inputSignal = SOFAspat(noiseSignal, hrtf,...
            azimuthAngleDeg(iAzimuth), elevationAngleDeg);
    
        % gammatone analysis filterbank - decompose signal into frequency bands
        subbandSignalArray = ...
            subbandDecompositionBinaural(inputSignal, FilterStates);
    
        % cumpute of interaural features
        [ipdRadCells, ivsCells, ~, ~] = ...
            interauralFeatureComputation(subbandSignalArray, ...
            AlgorithmParameters, FilterStates);
    
        ivsCellsCells{iAzimuth} = ivsCells;
    
        % calculate the mean over time of the IPD for each gammatone band
        for iBand = 1:nBands
            ipdRadMedianArray{iBand}(iAzimuth) = median(ipdRadCells{iBand});
        end
    end
    
    %% Fit the median IPD values to and create lookup struct
    lookuptable = cell(1, nBands);
    for iBand = 1:nBands
        [p, S, MU] = polyfit(ipdRadMedianArray{iBand}, azimuthAngleDeg,...
            polynomialOrder);
        lookuptable{iBand}.p = p;
        lookuptable{iBand}.MU = MU;
        lookuptable{iBand}.S = S;
    end
    
    %% for testing the IVS
    ivsArray = zeros(iAzimuth, iBand);
    for iAzimuth = 1:nAzimuth
        for iBand = 1:nBands
           ivsArray(iAzimuth, iBand) = ...
               numel(find(ivsCellsCells{iAzimuth}{iBand}));
        end
    end

end