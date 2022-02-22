% this is my on-the-fly adaptation of the dietz2011 code to fit my
% algorithm

%% Input signal
insig = repmat(sig_competingtalkers('five_speakers'),2,1);
% subsequently signal is split into frequency bands by auditory model in
% dietz2011.m, e.g.
% insig = auditory_model(insig);
fs = 16000;
flow = 50;
fhigh = 5000;

tau_cycles = 5;
%% Gammatone filtering
% find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(flow, fhigh, 1, []);

outsig = zeros(size(insig));
% gammatone filter centered at the center frequency for every frequency channel
for ii=1:length(fc)
  gammatone_filter = hohmann2002_filter (fs, fc(ii), ...
    fc(ii)/kv.fine_filter_finesse, ...
    kv.filter_attenuation_db, kv.filter_order);
  outsig(:,ii,1) = hohmann2002_process(gammatone_filter, squeeze(insig(:,1)));
  outsig(:,ii,2) = hohmann2002_process(gammatone_filter, squeeze(insig(:,2)));
end

%% IPD calculation
% ITF, eq. 2 in Dietz (2011)
outp.itf = outsig(:,:,2) .* conj(outsig(:,:,1));
% IPD, eq. 3 in Dietz (2011) but without low pass filtering
% outp.ipd = angle(outp.itf);
% LP filtering of IPD
tau = tau_cycles./fc;
a = exp( -1./(fs.*tau) );
% IPD of lowpassed signals, eq. 3 in Dietz (2011)
outp.ipd_lp = angle(lowpass(outp.itf,a));

% interaural coherence (IC) estimated by interaural vector strength (IVS), see
% eq. 7 in Dietz (2011)
outp.ic = interaural_vector_strength(outp.itf, tau, fs);

%% ILD LP filtering + calculation
% ILD filter
% low pass filter with a fixed cutoff frequency for every frequency channel
[b,a] = butter(2,30/(fs/2),'low');
outsig_ild = filter(b,a,outsig); % this uses insig in dietz2011 model

% interaural level difference, eq. 5 in Dietz (2011)
% max(sig,1e-4) avoids division by zero
outp.ild = 20.*log10(max(outsig_ild(:,:,2),1e-4)./max(outsig_ild(:,:,1),1e-4));

%% alternative version with only one type of LP filter (tau-based) for all
outsig_tau = lowpass(abs(outsig),a);
outsig_tau(abs(outsig_tau)<eps) = eps; % avoid division by zero and log(0)
outp.itf_tau = outsig_tau(:,:,2) .* conj(outsig_tau(:,:,1));
outp.ipd_tau = angle(outp.itf_tau);
outp.ild_tau = 20.*log10(outsig_tau(:,:,2)./outsig_tau(:,:,1));
% need seperate definition of IVS, too, where LP filtering is missing

%% ILD disambiguation of ITD (need to convert to write my own with IPD)
% convert interaural information into azimuth
itd_unwrapped = ...
    dietz2011_unwrapitd(outp.itd_lp,ild,outp.f_inst,2.5); % develop own

%% ITD to azimuth mapping (need to convert to IPD)
lookup = itd2angle_lookuptable(); % this is where my training signals go in
angl=itd2angle(itd_unwrapped,lookup);

%% this is itd2angle.m:
% phi = zeros(size(itd));
% for n = 1:size(itd,2)
%     % by calling the output S and MU, phi is z-scored, thus improving the fitting
%     phi(:,n)=polyval(lookup.p(:,n),itd(:,n),lookup.S{n},lookup.MU(:,n));
% end
% % neglect angles > 95°. Warning => maybe systematic underestimation for azi ~ 90°
% phi(abs(phi)>95) = NaN;

%% IVS decision mask application
angl_ivs_weighted = ...
    angl(outp.ic(:,n)>ic_threshold&[diff(outp.ic(:,n))>0; 0],n);

% histogram representation of detected angles
% h_ic=zeros(91,12);
% h_all=histc(angl,-90:2:90);
% for n=1:12
%     h_ic(:,n)=histc(angl(outp.ic(:,n)>ic_threshold&[diff(outp.ic(:,n))>0; 0],n),-90:2:90);
% end

%% ----- internal functions -----

% sig_competingtalkers
% auditoryfilterbank
% hohmann2002_filter and hohmann2002_process
% itd2angle_lookuptable
% itd2angle

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