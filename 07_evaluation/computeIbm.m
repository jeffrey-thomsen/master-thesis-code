% Compute ideal binary masks (IBM) for target and interferer glimpses by
% comparing signal levels of target and interferer signals
function [ibmTarget, ibmInterf] = ...
  computeIbm(targetSubbandSignals, interfSubbandSignals, AlgorithmParameters)

    % group subband signals into 23-ms time-frequency blocks and compute
    % energy

    energyTarget = zeros(size(targetSubbandSignals));
    energyInterf = zeros(size(targetSubbandSignals));

    winLen = round(0.023 * AlgorithmParameters.Gammatone.samplingRateHz);
    
    nWindows = floor(length(targetSubbandSignals)/winLen);
    nBands = AlgorithmParameters.Gammatone.nBands;
    for iWin = 1:nWindows
        tmpIndices = (1:winLen) + ((iWin-1)*winLen);
        for jBand = 1:nBands
            energyTarget(tmpIndices,jBand) = ...
                sum(abs(targetSubbandSignals(tmpIndices,jBand)).^2);
            energyInterf(tmpIndices,jBand) = ...
                sum(abs(interfSubbandSignals(tmpIndices,jBand)).^2);
        end
    end
    
    nRemainingSamples = mod(length(targetSubbandSignals),winLen)-1;
    for jBand = 1:nBands
        energyTarget(end-nRemainingSamples:end,jBand) = ...
            sum(abs(targetSubbandSignals(end-nRemainingSamples:end,jBand)).^2);
        energyInterf(end-nRemainingSamples:end,jBand) = ...
            sum(abs(interfSubbandSignals(end-nRemainingSamples:end,jBand)).^2);
    end

    % compute IBM from grouped energy time-frequency blocks
    snrCriterion = 0;
    SNR = 10 * log10(energyTarget ./ energyInterf);

    ibmTarget = SNR>snrCriterion;
    ibmInterf = SNR<-snrCriterion;

end