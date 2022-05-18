% plots a spectrogram from the mean of a stereo signal
% with parameters chosen for the DAGA publication
function [] = plotMonoSpectrogram(testSignal, samplingRateHz)

    meanTestSignal = (testSignal(:,1) + testSignal(:,2))./2;
    
    winlen = 2000;
    overlap = round(winlen-500);%kaiser5 *0.705);%blackmanharris *0.661);

    spectrogram(meanTestSignal, blackmanharris(winlen), overlap, [], ...
        samplingRateHz, 'yaxis', 'power');
    colormap('bone');
end