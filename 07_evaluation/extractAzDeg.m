% reshape azimuth estimate values into a matrix for easier data handling
function azDeg = extractAzDeg(ivsMask, azimuthDegCells, AlgorithmParameters, nBands)

    azDeg = zeros(nBands,length(ivsMask{1}));
    for iBand = 1:nBands
        if AlgorithmParameters.coherenceMask
            azDeg(iBand,ivsMask{iBand}) = azimuthDegCells{iBand};
        else
            azDeg(iBand,:) = azimuthDegCells{iBand};
        end
    end

end