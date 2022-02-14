function snr = calcSNRBinaural(sigma, delta)
        snr.L = calcSNR(sigma.L, delta.L);
        snr.R = calcSNR(sigma.R, delta.R);
end

function snr = calcSNR(sigma, delta)
    snr = sigma./(delta+eps);
end