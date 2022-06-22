% Calculate the CFR, defined as Sigma/Delta.
% Eps Offset added to avoid division by zero.
%
% Input:
% sigma, delta - structs containing Sigma and Delta signals of the left and
% right channel
%
% Output:
% cfr - struct containing the CFR values of the left and right channel
function cfr = calcCFRBinaural(sigma, delta)
        cfr.L = calcCFR(sigma.L, delta.L);
        cfr.R = calcCFR(sigma.R, delta.R);
end

function cfr = calcCFR(sigma, delta)
    cfr = sigma./(delta+eps);
end