function SNR = computeTimeFreqSnr(signalSubbandSignals, noiseSubbandSignals, ...
  samplingRateHz, nBands)

    % group into 23-ms time-frequency blocks and their compute energy
    energySignal = zeros(size(signalSubbandSignals));
    energyNoise = zeros(size(signalSubbandSignals));

    winLenSeconds = 0.023;
    winLenSamples = round(winLenSeconds * samplingRateHz);
    
    nWindows = floor(length(signalSubbandSignals)/winLenSamples);
    for iWin = 1:nWindows
        tmpIndices = (1:winLenSamples) + ((iWin-1)*winLenSamples);
        for jBand = 1:nBands
            energySignal(tmpIndices,jBand) = ...
                sum(abs(signalSubbandSignals(tmpIndices,jBand)).^2);
            energyNoise(tmpIndices,jBand) = ...
                sum(abs(noiseSubbandSignals(tmpIndices,jBand)).^2);
        end
    end
    
    nRemainingSamples = mod(length(signalSubbandSignals),winLenSamples)-1;
    for jBand = 1:nBands
        energySignal(end-nRemainingSamples:end,jBand) = ...
            sum(abs(signalSubbandSignals(end-nRemainingSamples:end,jBand)).^2);
        energyNoise(end-nRemainingSamples:end,jBand) = ...
            sum(abs(noiseSubbandSignals(end-nRemainingSamples:end,jBand)).^2);
    end

    % compute SNR from grouped energy time-frequency blocks
    SNR = 10 * log10(energySignal ./ energyNoise);



end