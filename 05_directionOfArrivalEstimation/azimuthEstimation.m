% Estimate the incidence angle of signal samples using an IPD- or 
% ITD-to-azimuth mapping function. Apply an IVS coherence mask to the
% azimuth estimates, passing only those whose interaural transfer function
% is sufficiently coherent.
%
% Input:
% ipdRadCells - cells containing IPD values of subband signals of all bands
% itdSecCells - cells containing ITD values of subband signals of all bands
% ivsMaskCells - cells containing logical arrays for each subband,
% determining which azimuth values to pass on for DOA estimation
% AlgorithmParameters - struct containing parametric information for the
% simulation, including the mapping polynomial coefficient lookup table
%
% Output:
% azimuthDegCells - cells for each subband containing the azimuth estimates
% in degrees that passed the IVS mask
function azimuthDegCells = azimuthEstimation(ipdRadCells, itdSecCells, ...
  ivsMaskCells, AlgorithmParameters)

    switch AlgorithmParameters.lookuptableType
        case 'itd'
            mappingData = itdSecCells;
        case 'ipd'
            mappingData = ipdRadCells;
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
            azimuthDeg = azimuthDeg(ivsMaskCells{iBand});
        elseif AlgorithmParameters.azimuthPooling
            if iBand>14
                azimuthDeg = median(cat(2, azimuthDegCells{1:14}), 2);
            end
        end

        % write computed values into structs for enhancement stage
        azimuthDegCells{iBand} = azimuthDeg;
    end
end