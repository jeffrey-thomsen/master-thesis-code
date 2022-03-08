% generate an IPD/ITD-to-azimuth mapping function in form of a polynomial
% hrtf - SOFA file
% 
% Input:
% hrtf - SOFA file containing binaural head-related impulse responses for
% various azimuth angles of incidence in the frontal hemisphere
% AlgorithmParameters - struct containing the algorithm parameters used for
% gammatone filterbank generation and interaural feature computation
% varargin - filename strings of sound files to be used for generating the
% mapping function. Multiple can be specified. White noise is used in case
% none are specified.
%
% Output:
% lookuptable - struct containing a set of polynomial coefficients for each
% gammatone band, to be used for mapping ITD or IPD values from the DOA
% estimation of the speech enhancement algorithm to azimuth values
% (estimated angle of incidence of a sound source)
function lookuptable = ...
  interauralToAzimuthLookuptable(hrtf, AlgorithmParameters, varargin)

    %% Configuration
    
    samplingRateHzHRTF = hrtf.Data.SamplingRate;
    AlgorithmParameters.Gammatone.samplingRateHz = samplingRateHzHRTF;

    % initialize filter states for interaural feature computation
    [FilterStates, nBands] = ...
        AlgorithmStatesConstructor(AlgorithmParameters);
    AlgorithmParameters.Gammatone.nBands = nBands;
    
    Obj = hrtf;

    polynomialOrder = 12; % Dietz2011 paper says 9;
    
    %% Preparation

    trainingSignal = [];
    if nargin > 2
        for iFile = 1:numel(varargin)
            % load training signal
            [signal, samplingRateHzSignal] = audioread(varargin{iFile}); %'sp01.wav');
            % adjust sampling rate to HRTF
            signal = resample(signal, samplingRateHzHRTF, samplingRateHzSignal);
            % append to training signal
            trainingSignal = [trainingSignal; signal];
        end
    else
        % length of noise signal used for the calculation (samples)
        nSamples = 1*samplingRateHzHRTF;
        % generate noise signal
        trainingSignal = randn(nSamples,1);
    end

    % normalize training signal
    trainingSignal = trainingSignal - mean(trainingSignal);
    trainingSignal = trainingSignal/max(abs(trainingSignal));
    
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

    azimuthAngleDeg = sort(azimuthAngleDeg);
    
    %% Calculation 
    
    nAzimuth = length(azimuthAngleDeg);
    
    ipdRadMedianArray = cell(1,nBands); %zeros(nAzimuth, nBands);
    itdSecMedianArray = cell(1,nBands); %zeros(nAzimuth, nBands);
%     ipdIvsMedianArray = cell(1,nBands); %zeros(nAzimuth, nBands);
    ivsCellsCells = cell(1, nAzimuth);
    
    for iAzimuth = 1:nAzimuth
    
        % generate noise coming from the given direction
        inputSignal = SOFAspat(trainingSignal, hrtf,...
            azimuthAngleDeg(iAzimuth), elevationAngleDeg);
    
        % gammatone analysis filterbank - decompose signal into frequency bands
        subbandSignalArray = ...
            subbandDecompositionBinaural(inputSignal, FilterStates);

        % compute interaural features
        [ipdRadCells, ivsCells, ~, ~, ~, itdSecCells, ~] = ...
            interauralFeatureComputation(subbandSignalArray, ...
            AlgorithmParameters, FilterStates);
    
        ivsCellsCells{iAzimuth} = ivsCells;
    
        % calculate the mean over time of the IPD/ITD for each gammatone band
        for iBand = 1:nBands
            ipdRadMedianArray{iBand}(iAzimuth) = median(ipdRadCells{iBand});
            itdSecMedianArray{iBand}(iAzimuth) = median(itdSecCells{iBand});
        end
    end
    
    %% Fit the median IPD values to and create lookup struct
    
    switch AlgorithmParameters.lookuptableType
        case 'itd'
            mappingData = itdSecMedianArray;
        case 'ipd'
            mappingData = ipdRadMedianArray;
    end

    lookuptable = cell(1, nBands);
    for iBand = 1:nBands
        [p, S, MU] = polyfit(mappingData{iBand}, azimuthAngleDeg,...
            polynomialOrder);
        lookuptable{iBand}.p = p;
        lookuptable{iBand}.MU = MU;
        lookuptable{iBand}.S = S;
    end
        
  
    %% for evaluation purposes

%     for i=[1,5,10,15,20]
%         figure(15); plot(azimuthAngleDeg,ipdRadMedianArray{i}); hold on;
% %         figure; plot(azimuthAngleDeg,ipdRadDisMedianArray{i})
% %         hold on; plot(azimuthAngleDeg,unwrap(ipdRadMedianArray{i}));
% %         plot(azimuthAngleDeg,(ipdRadMedianArray{i}));
% %         plot(azimuthAngleDeg,unwrap(ipdRadMedianArray{i})-2*pi);
%     end
% 
%     for iBand=1:30
%         figure(16); 
%         plot(ipdRadMedianArray{iBand},azimuthAngleDeg); 
%         hold on; 
%         ylim([-100, 100]);
%         xlim([-pi, pi]);
%         figure(17); 
%         plot(linspace(-pi, pi), polyval(lookuptable{iBand}.p,...
%             linspace(-pi, pi), lookuptable{iBand}.S, lookuptable{iBand}.MU)); 
%         hold on; 
%         ylim([-100, 100]);
%         xlim([-pi, pi]);
%     end
% % 
% %     for iBand=[1,5,10,15,20]
% %         figure(18); plot(ildMedianArray{iBand},azimuthAngleDeg); hold on; ylim([-100, 100]); xlim([-10, 10]);
% %         figure(19); plot(linspace(-10,10),polyval(pILD{iBand}, linspace(-pi,pi), SILD{iBand}, ...
% %             MUILD{iBand})); hold on; ylim([-100, 100]); xlim([-10, 10]);
% %     end

    %% for testing the IVS
%     ivsArray = zeros(iAzimuth, iBand);
%     for iAzimuth = 1:nAzimuth
%         for iBand = 1:nBands
%            ivsArray(iAzimuth, iBand) = ...
%                numel(find(ivsCellsCells{iAzimuth}{iBand}));
%         end
%     end

end