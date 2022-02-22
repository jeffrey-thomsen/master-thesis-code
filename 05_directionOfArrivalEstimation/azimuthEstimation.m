function azimuthDegCells = azimuthEstimation(ipdRad, ...
  ivsMask, AlgorithmParameters)

    nBands = AlgorithmParameters.Gammatone.nBands;
    azimuthDegCells = cell(1, nBands);

    for iBand = 1:nBands

        % evaluate IPD-to-azimuth mapping polynomial
        azimuthDeg = ipdToAzimuthMapping(ipdRad{iBand}, ...
            AlgorithmParameters.lookuptable{iBand});
        
        % apply binary IVS filter mask - remove all azimuth estimates for 
        % which coherence criteria were not met
        azimuthDeg = azimuthDeg(ivsMask{iBand});
        
        % write computed values into structs for enhancement stage
        azimuthDegCells{iBand} = azimuthDeg;
    end
end