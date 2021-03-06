% Calculate the interaural coherence (IC) with the help of the interaural vector
% strength (IVS) as defined by eq. 7 in Dietz (2011)
% author: Mathias Dietz, adapted by Jeffrey Thomsen
%
% Input:
% itf - Nx1 complex-valued interaural transfer function vector
% tauS - time constant for filtering
% samplingRateHz - sampling rate of the signal in Hz
% ivsThreshold - threshold value for binary IVS mask
% States - filter states to be handed to the next signal block processed
%
% Output:
% ivsMask - logical array specifiyng which signal samples passed the
% coherence mask (threshold and rising slope conditions)
% States - see above
function [ivsMask, States] = calcIvs(itf, tauS, samplingRateHz, ...
  ivsThreshold, States)

    a = exp(-1./(samplingRateHz.*tauS));

    [numerator, States.ivsNumeratorLP] = ...
        filter(1-a, [1 -a], itf, States.ivsNumeratorLP);
    [denominator, States.ivsDenominatorLP] = ...
        filter(1-a, [1 -a], abs(itf), States.ivsDenominatorLP);

    ivs = abs(numerator)./abs(denominator);
    
    % threshold and rising slope conditions
    ivsMask = ...
        (ivs >= ivsThreshold) & ...
        (diff([States.ivsPreviousValue; ivs]) >= 0);

    States.ivsPreviousValue = ivs(end);

end