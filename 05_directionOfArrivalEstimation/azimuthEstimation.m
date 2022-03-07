function azimuthDegCells = azimuthEstimation(ipdRad, itdSec, ...
  ivsMask, AlgorithmParameters)

    switch AlgorithmParameters.lookuptableType
        case 'itd'
            mappingData = itdSec;
        case 'ipd'
            mappingData = ipdRad;
    end
    
    if AlgorithmParameters.coherenceMask && AlgorithmParameters.azimuthPooling
        warning('coherenceMask and azimuthPooling both set to True, will only apply coherenceMask')
    end

    nBands = AlgorithmParameters.Gammatone.nBands;
    azimuthDegCells = cell(1, nBands);
    
    for iBand = 1:nBands

        % evaluate IPD/ITD-to-azimuth mapping polynomial

        azimuthDeg = interauralToAzimuthMapping(mappingData{iBand}, ...
            AlgorithmParameters.lookuptable{iBand});

        if AlgorithmParameters.coherenceMask
            % apply binary IVS filter mask - remove all azimuth estimates for 
            % which coherence criteria were not met
            azimuthDeg = azimuthDeg(ivsMask{iBand});
        elseif AlgorithmParameters.azimuthPooling
            if iBand>14
                azimuthDeg = median(cat(2,azimuthDegCells{1:14}),2);
            end
        end

        % write computed values into structs for enhancement stage
        azimuthDegCells{iBand} = azimuthDeg;
    end
end