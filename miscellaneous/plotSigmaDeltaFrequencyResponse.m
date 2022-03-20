% for fs=1000 Hz and a period of k=5 samples, equiv. to 0.005 s or 200 Hz,
% plot the frequency responses of
% Sigma = x[n] + x[n-5]
% Delta = x[n] - x[n-5]

figure;
hold on;

[H,W] = freqz([1 0 0 0 0 +1], 1, 10000, 1000);
plot(W,10*log10(abs(H)));

[H,W] = freqz([1 0 0 0 0 -1], 1, 10000, 1000);
plot(W,10*log10(abs(H)));

xlim([0 500])
ylim([-20 5])
legend('\Sigma enhancement', '\Delta cancellation')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
set(gca,'FontSize',18);