AlgorithmParameters = AlgorithmParametersConstructor();
[AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
    AlgorithmStatesConstructor(AlgorithmParameters);

fs = AlgorithmStates.L.GammatoneStates.analyzer.fs;

freqvector = 1:0.1:10000;
zfreqvector = exp(2*1i*pi*freqvector/fs);
freqresponse = hohmann2002_freqz(AlgorithmStates.L.GammatoneStates.analyzer, zfreqvector);

figure('Position',[100 100 2*576 432])
for iBand = 1:2:30
semilogx(freqvector,20*log10(abs(freqresponse(:,iBand)))-6,'b','LineWidth',2)
hold on
end
xlim([10 8000])
ylim([-40 1])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
set(gca,'FontSize',16);