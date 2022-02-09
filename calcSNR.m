function snr = calcSNR(sigma, delta)
    snr = sigma./(delta+eps);
end