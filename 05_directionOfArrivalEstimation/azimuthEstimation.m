function azimuthDegCells = azimuthEstimation(ipdRad, itdSec, ...
  ivsMask, AlgorithmParameters)

    switch AlgorithmParameters.lookuptableType
        case 'itd'
            mappingData = itdSec;
        case 'ipd'
            mappingData = ipdRad;
    end
    
    nBands = AlgorithmParameters.Gammatone.nBands;
    azimuthDegCells = cell(1, nBands);
    
    for iBand = 1:nBands

        % evaluate IPD-to-azimuth mapping polynomial

        azimuthDeg = interauralToAzimuthMapping(mappingData{iBand}, ...
            AlgorithmParameters.lookuptable{iBand});

        % evaluate ITD-to-azimuth mapping polynomial
%         azimuthDeg = ipdToAzimuthMapping(itdSec{iBand}, ...
%             AlgorithmParameters.itdLookuptable{iBand});
        
        if AlgorithmParameters.coherenceMask
            % apply binary IVS filter mask - remove all azimuth estimates for 
            % which coherence criteria were not met
            azimuthDeg = azimuthDeg(ivsMask{iBand});
        end
        
        % write computed values into structs for enhancement stage
        azimuthDegCells{iBand} = azimuthDeg;
    end
end