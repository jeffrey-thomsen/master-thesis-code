% Compute ideal binary masks (IBM) for target and interferer glimpses by
% comparing signal levels of target and interferer signals
function [ibmTarget, ibmInterf] = ...
  computeIbm(targetSubbandSignals, interfSubbandSignals, AlgorithmParameters)
    
    nBands = AlgorithmParameters.Gammatone.nBands;
    
    % lowpass filter subband signals for some temporal smoothing
    targetSubbandSignals = firstOrderLowPass(targetSubbandSignals, ...
            zeros(1, nBands), AlgorithmParameters);
    interfSubbandSignals = firstOrderLowPass(interfSubbandSignals, ...
            zeros(1, nBands), AlgorithmParameters);
    
    % compute IBM from LP-filtered subband signals
    snrCriterion = 0;
    SNR = 10 * log10(abs(targetSubbandSignals).^2 ./ abs(interfSubbandSignals).^2);

    ibmTarget = SNR>snrCriterion;
    ibmInterf = SNR<-snrCriterion;

end