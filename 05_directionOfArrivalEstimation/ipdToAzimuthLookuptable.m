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
    [noiseSignal,~]=audioread('sp01.wav');
    [noiseSignal2,~]=audioread('sp30.wav');
    noiseSignal = [noiseSignal; noiseSignal2];
    
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
    ipdRadDisMedianArray = cell(1,nBands); %zeros(nAzimuth, nBands);
%     ipdIvsMedianArray = cell(1,nBands); %zeros(nAzimuth, nBands);
    ivsCellsCells = cell(1, nAzimuth);
    
    for iAzimuth = 1:nAzimuth
    
        % generate noise coming from the given direction
        inputSignal = SOFAspat(noiseSignal, hrtf,...
            azimuthAngleDeg(iAzimuth), elevationAngleDeg);
    
        % gammatone analysis filterbank - decompose signal into frequency bands
        subbandSignalArray = ...
            subbandDecompositionBinaural(inputSignal, FilterStates);
    
        % compute of interaural features
        [ipdRadCells, ivsCells, ildCells, ~, ~, itdCells, ~] = ...
            interauralFeatureComputation(subbandSignalArray, ...
            AlgorithmParameters, FilterStates);
    
        ivsCellsCells{iAzimuth} = ivsCells;
    
        % calculate the mean over time of the IPD for each gammatone band
        for iBand = 1:nBands
            ipdRadMedianArray{iBand}(iAzimuth) = median(ipdRadCells{iBand});
            ipdRadDis = disambiguateIpd(ipdRadCells{iBand}, ildCells{iBand});
            ipdRadDisMedianArray{iBand}(iAzimuth) = median(ipdRadDis);
            ildMedianArray{iBand}(iAzimuth) = median(ildCells{iBand});
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

%         [pILD{iBand}, SILD{iBand}, MUILD{iBand}] = polyfit(ildMedianArray{iBand}, azimuthAngleDeg,...
%             polynomialOrder);
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
%         figure(16); plot(ipdRadMedianArray{iBand},azimuthAngleDeg); hold on; ylim([-100, 100]); xlim([-pi, pi]);
%         figure(17); plot(linspace(-pi,pi),polyval(lookuptable{iBand}.p, linspace(-pi,pi), lookuptable{iBand}.S, ...
%             lookuptable{iBand}.MU)); hold on; ylim([-100, 100]); xlim([-pi, pi]);
%     end
% 
% %     for iBand=[1,5,10,15,20]
% %         figure(18); plot(ildMedianArray{iBand},azimuthAngleDeg); hold on; ylim([-100, 100]); xlim([-10, 10]);
% %         figure(19); plot(linspace(-10,10),polyval(pILD{iBand}, linspace(-pi,pi), SILD{iBand}, ...
% %             MUILD{iBand})); hold on; ylim([-100, 100]); xlim([-10, 10]);
% %     end

    %% for testing the IVS
    ivsArray = zeros(iAzimuth, iBand);
    for iAzimuth = 1:nAzimuth
        for iBand = 1:nBands
           ivsArray(iAzimuth, iBand) = ...
               numel(find(ivsCellsCells{iAzimuth}{iBand}));
        end
    end

end