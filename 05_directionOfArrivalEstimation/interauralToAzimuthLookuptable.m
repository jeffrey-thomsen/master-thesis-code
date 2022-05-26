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
        hrtf.SourcePosition(:,2) == elevationAngleDeg & ... % horizontal plane
        (((hrtf.SourcePosition(:,1) >= -90 & hrtf.SourcePosition(:,1) < 0) ... % front right
        | hrtf.SourcePosition(:,1) >= 270) | ... % alternative front right
        (hrtf.SourcePosition(:,1) >= 0 & hrtf.SourcePosition(:,1) <= 90)); % front left
    azimuthAngleDeg = hrtf.SourcePosition(idFrontHor,1);
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
%     centerFreqsHz = FilterStates.L.GammatoneStates.analyzer.center_frequencies_hz;
% 
%     switch AlgorithmParameters.lookuptableType
%         case 'itd'
%             mappingSpace = linspace(-1e-3, 1e-3, 1000);
%         case 'ipd'
%             mappingSpace = linspace(-pi, pi, 1000);
%     end
% 
%     i = 1;
%     for iBand=1:4:17
%         figure(1);
% %         set(gca,'FontSize',18);
%         plot(mappingData{iBand}, azimuthAngleDeg); 
%         hold on;
%         ylim([-90, 90]);
%         title('Data for fitting ITD-to-azimuth polynomial')
%         xlabel('ITD (s)')
%         ylabel('azimuth (degrees)')
% 
%         figure(2);
% %         set(gca,'FontSize',18);
%         plot(mappingSpace, polyval(lookuptable{iBand}.p, mappingSpace, ...
%             lookuptable{iBand}.S, lookuptable{iBand}.MU)); 
%         hold on;
%         ylim([-90, 90]);
%         title('Fitted ITD-to-azimuth polynomial')
%         xlabel('ITD (s)')
%         ylabel('azimuth (degrees)')
% 
%         legendString{i}=['fc=',num2str(centerFreqsHz(iBand),'%.0f'),'Hz'];
%         i = i+1;
%     end
%     figure(1); legend(legendString,'Location','southeast')
%     figure(2); legend(legendString,'Location','southeast')
% 
%         i = 1;
%     for iBand=18:4:30
%         figure(3);
% %         set(gca,'FontSize',18);
%         plot(mappingData{iBand}, azimuthAngleDeg); 
%         hold on;
%         ylim([-90, 90]);
%         title('Data for fitting ITD-to-azimuth polynomial')
%         xlabel('ITD (s)')
%         ylabel('azimuth (degrees)')
% 
%         figure(4); 
% %         set(gca,'FontSize',18);
%         plot(mappingSpace, polyval(lookuptable{iBand}.p, mappingSpace, ...
%             lookuptable{iBand}.S, lookuptable{iBand}.MU)); 
%         hold on;
%         ylim([-90, 90]);
%         title('Fitted ITD-to-azimuth polynomial')
%         xlabel('ITD (s)')
%         ylabel('azimuth (degrees)')
% 
%         legendString{i}=['fc=',num2str(centerFreqsHz(iBand),'%.0f'),'Hz'];
%         i = i+1;
%     end
%     figure(3); legend(legendString,'Location','southeast')
%     figure(4); legend(legendString,'Location','southeast')


    %% for testing the IVS
%     ivsArray = zeros(iAzimuth, iBand);
%     for iAzimuth = 1:nAzimuth
%         for iBand = 1:nBands
%            ivsArray(iAzimuth, iBand) = ...
%                numel(find(ivsCellsCells{iAzimuth}{iBand}));
%         end
%     end

end