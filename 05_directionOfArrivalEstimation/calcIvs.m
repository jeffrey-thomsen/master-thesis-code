% Calculate the interaural coherence (IC) with the help of the interaural vector
% strength (IVS) as defined by eq. 7 in Dietz (2011)
function [ivsMask, States] = calcIvs(itf,tauS,samplingRateHz,ivsThreshold,States)

    % tau_coherence = 15e-3; % good value for ipd_fine
    a = exp(-1./(samplingRateHz.*tauS));

    [numerator, States.ivsNumeratorLP] = ...
        filter(1-a, [1 -a], itf, States.ivsNumeratorLP);
    [denominator, States.ivsDenominatorLP] = ...
        filter(1-a, [1 -a], abs(itf), States.ivsDenominatorLP);
    %     when is this identical to my simple lowpass?? tau=1/fc or
    %     1/(2*pi*fc)???
    %     [numerator, States.ivsNumeratorLP] = ...
    %         firstOrderLowPass(itf, States.ivsNumeratorLP ,AlgorithmParameters, 1/30);
    %     [denominator, States.ivsDenominatorLP] = ...
    %         firstOrderLowPass(abs(itf), States.ivsDenominatorLP ,AlgorithmParameters, 1/30);

    ivs = abs(numerator)./abs(denominator);
    % threshold and rising slope conditions
    ivsMask = ...
        (ivs >= ivsThreshold) & ...
        (diff([States.ivsPreviousValue; ivs]) >= 0);

    States.ivsPreviousValue = ivs(end);

end