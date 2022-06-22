% plot heat map of time-frequency TIR improvement values
function [] = plotTimeFreqBlockTirHeatmap(TIR, samplingRateHz, centerFreqsHz) 
    
    winLenSeconds = 0.023;
    winLenSamples = round(winLenSeconds * samplingRateHz);
    nWindows = ceil(length(TIR)/winLenSamples);
    
    figure;
    heatmap((1:nWindows)*winLenSeconds, flip(centerFreqsHz), ...
        flip(TIR(1:winLenSamples:end,:)'))
    annotation('textarrow', [0.97,0.97], [0.5,0.5], 'string', 'SNR (dB)', ...
      'HeadStyle', 'none', 'LineStyle', 'none', 'HorizontalAlignment', ...
      'center', 'TextRotation', 90);
    title('SNR heatmap');
    xlabel('Time (s)')
    ylabel('Gammatone band center frequency (Hz)')

end