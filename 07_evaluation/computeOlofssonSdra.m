% Hilbert-transform based measure of nonlinear distortion in speech
% enhancement algorithm, expressed as A-weighted signal-to-distortion ratio
% (SDR_A), according to Olofsson2006
function psd_delta_total_A_dB = computeOlofssonSdra(testSignal, ...
  targetAngle, AlgorithmParameters, AlgorithmStates , Plot)
    
    % compute Hilbert transform
    testSignal_analytical = hilbert(testSignal);
    testSignal_hilbert = imag(testSignal_analytical);
    
    % process test signal and its Hilbert transform separately
    AlgorithmParameters.targetRangeDeg = targetAngle + [-5 5];

    [enhancedSignal, ~, ~] = ...
        speechEnhancement(testSignal, ...
        AlgorithmParameters, AlgorithmStates);

    [enhancedSignal_hilb, ~, ~] = ...
        speechEnhancement(testSignal_hilbert, ...
        AlgorithmParameters, AlgorithmStates);
    
    % construct enhanced analytical signal
    enhancedSignal_analytical = enhancedSignal + 1j*enhancedSignal_hilb;
    
    % compute power spectral density (PSD)
    noverlap = 512;
    [enhanced_psd, f] = pwelch(enhancedSignal_analytical, ...
        hamming(2*noverlap), noverlap, [], ...
        AlgorithmParameters.Gammatone.samplingRateHz);
    enhanced_psd_dB = 10*log10(enhanced_psd);
    
    % decompose PSD into positive and negative frequencies and
    psd_signal_dB = enhanced_psd_dB(1:noverlap,:);
    psd_distortion_dB = enhanced_psd_dB(end:-1:noverlap+1,:);
    
    if Plot
        figure;
        semilogx(f(1:noverlap),psd_signal_dB);
        hold on;
        semilogx(f(1:noverlap),psd_distortion_dB);
    end
    
    %
    psd_delta_dB = psd_signal_dB - psd_distortion_dB;
    psd_delta_total_dB = 10*log10((1/numel(psd_delta_dB))*...
                                sum(10.^(0.1*psd_delta_dB), "all"));
    
    % A Weighting curve
    RA = (12194^2 * f.^4) ./ ((f.^2 + 20.6^2) .* ...
        sqrt((f.^2+107.7^2).*(f.^2+737.9^2)) .* (f.^2 + 12194^2));
    
    AWeighting_dB = 20*log10(RA)+2;
    
    psd_delta_A_dB = psd_delta_dB + AWeighting_dB(1:noverlap);
    
    psd_delta_total_A_dB = 10*log10((1/numel(psd_delta_A_dB))*...
                                sum(10.^(0.1*psd_delta_A_dB), "all"));

end