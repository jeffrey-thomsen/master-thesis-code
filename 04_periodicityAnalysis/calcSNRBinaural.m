% Calculate the SNR, defined as Sigma/Delta.
% Eps Offset added to avoid division by zero.
%
% Input:
% sigma, delta - structs containing Sigma and Delta signals of the left and
% right channel
%
% Output:
% snr - struct containing the SNR values of the left and right channel
function snr = calcSNRBinaural(sigma, delta)
        snr.L = calcSNR(sigma.L, delta.L);
        snr.R = calcSNR(sigma.R, delta.R);
end

function snr = calcSNR(sigma, delta)
    snr = sigma./(delta+eps);
end