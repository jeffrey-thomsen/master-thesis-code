%% from exp_dietz2011.m

signal=repmat(sig_competingtalkers('five_speakers'),2,1);
fs = 16000;
signal=signal(1:5*fs,:);
s_pos = [-80 -30 0 30 80];
ic_threshold=0.98;
panellabel = 'ab';

% run IPD model on signal
[hairc_fine, fc, hairc_ild]=dietz2011(signal,fs,'fhigh',1400,flags.disp);

% convert interaural information into azimuth
itd_unwrapped = ...
    dietz2011_unwrapitd(hairc_fine.itd_lp,hairc_ild,hairc_fine.f_inst,2.5);
angl=itd2angle(itd_unwrapped,lookup);

h_ic=zeros(91,12);
h_all=histc(angl,-90:2:90);
for n=1:12
    h_ic(:,n)=histc(angl(hairc_fine.ic(:,n)>ic_threshold&[diff(hairc_fine.ic(:,n))>0; 0],n),-90:2:90);
end
example_output.angle_fine = angl;
example_output.IVS_fine = hairc_fine.ic;
example_output.histogram_angle_label = -90:2:90;
example_output.histograms_with_IVS = h_ic;
example_output.histograms_without_IVS = h_all;

%% from dietz2011.m

function [fine,fc,ild,env] = dietz2011(insig,fs,varargin)
    % ---- modulation filterbank ------
    % filter signals with three different filters to get fine structure, envelope
    % and ILD low pass
    amt_disp('  apply second filterbank',flags.disp)
    [inoutsig_fine,fc_fine,inoutsig_env,fc_env,inoutsig_ild] = ...
      dietz2011_filterbank(inoutsig,fs,fc,'argimport',flags,kv);
    
    % ---- binaural processor ------
    % calculate interaural parameters for fine structure and envelope and calculate
    % ILD
    % -- fine structure
    amt_disp('  calculating interaural functions from haircell fine structure',flags.disp);
    fine = dietz2011_interauralfunctions(inoutsig_fine,fs,fc_fine,'argimport',flags,kv);
    % --envelope
    amt_disp('  calculating interaural functions from haircell modulation',flags.disp);
    env = dietz2011_interauralfunctions(inoutsig_env,fs, ...
      kv.mod_center_frequency_hz+0*fc_env,'argimport',flags,kv);
    % -- ILD
    % interaural level difference, eq. 5 in Dietz (2011)
    % max(sig,1e-4) avoids division by zero
    amt_disp('  determining ILD',flags.disp);
    ild = 20/kv.compression_power*log10(max(inoutsig_ild(:,:,2),1e-4)./max(inoutsig_ild(:,:,1),1e-4));
end

%% arg_dietz2011_interauralfunctions.m

function definput=arg_dietz2011_interauralfunctions(definput)

  definput.keyvals.signal_level_dB_SPL = 70;
  definput.keyvals.tau_cycles  = 5; % see Fig. 3c in Dietz (2011)
  definput.keyvals.compression_power = 0.4;

  % ask for simulating temporal resolution of binaural processor
  % this returns the *_lp values
  definput.flags.lowpass = {'lowpass','nolowpass'};
end

%% from dietz2011_interauralfuncions.m

function [outp] = dietz2011_interauralfunctions(insig,fs,fc,varargin)

definput.import = {'dietz2011_interauralfunctions'};
[flags,kv]  = ltfatarghelper({},definput,varargin);

% ----- interaural parameters -----
% interaural transfer function (ITF), eq. 2 in Dietz (2011)
outp.itf = insig(:,:,2) .* conj(insig(:,:,1));
% interaural phase difference (IPD), eq. 3 in Dietz (2011) but without low pass
% filtering
outp.ipd = angle(outp.itf);
% ----- low pass versions of interaural parameters -----
% The low pass is simulating a finite time resolution of the binaural system and
% is given by the coh_cycles parameter
tau = kv.tau_cycles./fc;
if flags.do_lowpass
  a = exp( -1./(fs.*tau) );
  % interaural phase difference (IPD) of lowpassed signals, eq. 3 in Dietz (2011)
  outp.ipd_lp = angle(lowpass(outp.itf,a));
  outp.f_inst_lp = lowpass(outp.f_inst,a);
  outp.itd_lp = 1/(2*pi)*outp.ipd_lp./outp.f_inst_lp;
  for k = 1:length(fc)
    outp.itd_C_lp(:,k) = 1/(2*pi)*outp.ipd_lp(:,k)/fc(k);
  end
  % interaural level difference
  inoutsig = lowpass(abs(insig),a);
  inoutsig(abs(inoutsig)<eps) = eps; % avoid division by zero and log(0)
  outp.ild_lp = 20./kv.compression_power.*log10(inoutsig(:,:,2)./inoutsig(:,:,1));
end

% interaural coherence (IC) estimated by interaural vector strength (IVS), see
% eq. 7 in Dietz (2011)
outp.ic = interaural_vector_strength(outp.itf, tau, fs);

% weighting of channels for cumulative ixd determination
% sqrt(2) is due to half-wave rectification (included 28th Sep 07)
outp.rms = kv.signal_level_dB_SPL*kv.compression_power + 20*log10(sqrt(2)*min(rms(squeeze(abs(insig(:,:,1)))),rms(squeeze(abs(insig(:,:,2))))));
outp.rms = max(outp.rms,eps); % avoid negative weights

end

% ----- internal functions -----
% lowpass
function outsig = lowpass(sig, a)
  % outsig = lowpass(sig, a)
  % This is a simple low-pass filter y(n) = (1-a)*x(n) - a*y(n-1)
  % Meaning of parameter a:
  %   a - damping coefficient (0 - no filtering, 1 - flat output)
  %   tau = 1/(2*pi*f_c)      where f_c is the cutoff frequency of the filter
  %   a = exp(-1/(fs*tau))   where fs - sampling frequency
  %
  if ndims(sig)==3
    sig = permute(sig,[1 3 2]); % => [samples ears channels]
    channels = size(sig,3);
    outsig = zeros(size(sig));
    for ii=1:channels
      outsig(:,:,ii) = filter([1-a(ii)], [1, -a(ii)], sig(:,:,ii));
    end
    outsig = permute(outsig,[1 3 2]); % => [samples channels ears]
  else
    channels = size(sig,2);
    outsig = zeros(size(sig));
    for ii=1:channels
      outsig(:,ii) = filter([1-a(ii)], [1, -a(ii)], sig(:,ii));
    end
  end
end

% Calculate the interaural coherence (IC) with the help of the interaural vector
% strength (IVS) as defined by eq. 7 in Dietz (2011)
function ivs = interaural_vector_strength(itf,tau_coherence,fs)
  %tau_coherence = 15e-3; % good value for ipd_fine
  c_coh = exp(-1./(fs.*tau_coherence));
  if length(tau_coherence)==1
    ivs = abs(filter(1-c_coh,[1 -c_coh],itf))./abs(filter(1-c_coh,[1 -c_coh],abs(itf)));
  elseif length(tau_coherence)==size(itf,2)
    ivs = zeros(size(itf));
    for n=1:length(tau_coherence)
      ivs(:,n) = abs(filter(1-c_coh(n),[1 -c_coh(n)],itf(:,n)))./ ...
          abs(filter(1-c_coh(n),[1 -c_coh(n)],abs(itf(:,n))));
    end
  else
    error('wrong number of tau_coherence values')
  end
end

%% from dietz2011_filterbank.m

% Split the input signal:
% - below 1.4 kHz the fine structure filter is applied
% - above 1.4 kHz the modulation filter is applied
fc_fine = fc(fc<=1400);
fc_env = fc(fc>1400);
insig_fine = insig(:,fc<=1400,:);
insig_env = insig(:,fc>1400,:);

% --- fine structur filter ---
outsig_fine = zeros(size(insig_fine));
% gammatone filter centered at the center frequency for every frequency channel
for ii=1:length(fc_fine)
  gammatone_filter = hohmann2002_filter (fs, fc_fine(ii), ...
    fc_fine(ii)/kv.fine_filter_finesse, ...
    kv.filter_attenuation_db, kv.filter_order);
  outsig_fine(:,ii,1) = hohmann2002_process(gammatone_filter, squeeze(insig_fine(:,ii,1)));
  outsig_fine(:,ii,2) = hohmann2002_process(gammatone_filter, squeeze(insig_fine(:,ii,2)));
end

% --- modulation/envelope filter ---
outsig_env = zeros(size(insig_env));
% gammatone filter centered at a fixed frequency for every frequency channel
gammatone_filter = hohmann2002_filter (fs, kv.mod_center_frequency_hz, ...
  kv.mod_center_frequency_hz/kv.mod_filter_finesse, ...
  kv.filter_attenuation_db, kv.filter_order);
for ii=1:length(fc_env)
  outsig_env(:,ii,1) = hohmann2002_process(gammatone_filter, squeeze(insig_env(:,ii,1)));
  outsig_env(:,ii,2) = hohmann2002_process(gammatone_filter, squeeze(insig_env(:,ii,2)));
end

% --- ILD filter ---
% low pass filter with a fixed cutoff frequency for every frequency channel
[b,a] = butter(kv.level_filter_order,kv.level_filter_cutoff_hz/(fs/2),'low');
outsig_ild = filter(b,a,insig);