% This script has been adapted from the Auditory Modeling Toolbox to demo
% the Hohmann2002 gammatone filterbank working only with the code used in 
% the Thomsen2022 speech enhancement algorithm. The original can bo found
% under the following link:
% http://amtoolbox.org/amt-1.1.0/doc/demos/demo_hohmann2002.php
%
% Reference:
% Piotr Majdak, Clara Hollomey, and Robert Baumgartner. "AMT 1.0: The 
% toolbox for reproducible research in auditory modeling." submitted to 
% Acta Acustica (2021).
% https://amtoolbox.org/

%%
% Shows how to use the gammatone filterbank from Hohmann(2002)
%
%   This Example demonstrates how to create and how to use the
%   combined analysis-synthesis Filterbank system.
%
%   Figure 4: Figure 4 shows the impulse response of the analysis-synthesis system in the time domain
%
%   Figure 5: Figure 5 shows shows its frequency response
%
%   See also: exp_hohmann2002 hohmann2002 hohmann2002_process


%% First, create a filterbank analyzer %%%

flow =    70;
fhigh =  6700;
base_frequency_hz         =  1000;
sampling_rate_hz          = 16000;
filters_per_ERB           =     1.0;
desired_delay_in_seconds  =     0.004;
filter_order              =     4;
bandwidth_factor          =     1.0;

disp('Building analysis filterbank');
analyzer = hohmann2002(sampling_rate_hz, flow, base_frequency_hz, fhigh, ...
    filters_per_ERB, filter_order, bandwidth_factor);


%% Now create a synthesizer that can resynthesize the analyzer's output %%%

disp(['Building synthesizer for an analysis-synthesis delay of ', num2str(desired_delay_in_seconds), ' seconds']);
synthesizer = hohmann2002_synth(analyzer, desired_delay_in_seconds);

%% Extract the synthesizer's parameters %%%
disp(['The synthesizers parameters:','----------------------------']);
delay = synthesizer.delay;
mixer = synthesizer.mixer;

bands = length(mixer.gains);

disp(sprintf('%3s|%7s | %22s | %5s', '# ', 'delay ', 'phase factor    ', 'gain / dB'));

for band = 1:bands
  disp(fprintf('%3d|%7d | %9f + %9fi | %5.2f', band,delay.delays_samples(band),real(delay.phase_factors(band)), imag(delay.phase_factors(band)), 20*log10(mixer.gains(band))));
end

%%  plot the resynthesized impulse and the frequency response of the  %%%
%%  analysis-synthesis system                                         %%%

impulse = [1; zeros(8191,1)];
[analyzed_impulse, analyzer] = hohmann2002_process(analyzer, impulse);
[resynthesized_impulse, synthesizer] = hohmann2002_process(synthesizer, analyzed_impulse);

figure(4);
plot([0:8191]/sampling_rate_hz*1e3, resynthesized_impulse);
axis([40/sampling_rate_hz*1e3, 120/sampling_rate_hz*1e3, -1, 1]);
title('impulse response of the analysis-synthesis system');
xlabel('time / ms');
ylabel('system output');


disp('Figure 4 shows the impulse response of the analysis-synthesis');
disp('system in the time domain.');
disp('Figure 5 shows its frequency response.');

frequency = [0:8191] * sampling_rate_hz / 8192;
figure(5)
plot(frequency, 20 * log10(abs(fft(resynthesized_impulse'))));
axis([0, sampling_rate_hz/2, -40, 5]);
title('frequency response of the analysis-synthesis-system');
xlabel('frequency / Hz');
ylabel('system response level / dB'); 