% compute TIR in time-frequency blocks for generation of an IBM
% target and interferer subband signals are grouped into time-frequency
% blocks and the TIR is computed for each block
function TIR = computeTimeFreqTir(targetSubbandSignals, interfSubbandSignals, ...
  samplingRateHz, nBands)

    % group into 23-ms time-frequency blocks and their compute energy
    energyTarget = zeros(size(targetSubbandSignals));
    energyInterf = zeros(size(targetSubbandSignals));

    winLenSeconds = 0.023;
    winLenSamples = round(winLenSeconds * samplingRateHz);
    
    nWindows = floor(length(targetSubbandSignals)/winLenSamples);
    for iWin = 1:nWindows
        tmpIndices = (1:winLenSamples) + ((iWin-1)*winLenSamples);
        for jBand = 1:nBands
            energyTarget(tmpIndices,jBand) = ...
                sum(abs(targetSubbandSignals(tmpIndices,jBand)).^2);
            energyInterf(tmpIndices,jBand) = ...
                sum(abs(interfSubbandSignals(tmpIndices,jBand)).^2);
        end
    end
    
    nRemainingSamples = mod(length(targetSubbandSignals),winLenSamples)-1;
    for jBand = 1:nBands
        energyTarget(end-nRemainingSamples:end,jBand) = ...
            sum(abs(targetSubbandSignals(end-nRemainingSamples:end,jBand)).^2);
        energyInterf(end-nRemainingSamples:end,jBand) = ...
            sum(abs(interfSubbandSignals(end-nRemainingSamples:end,jBand)).^2);
    end

    % compute TIR from grouped energy time-frequency blocks
    TIR = 10 * log10(energyTarget ./ energyInterf);



end