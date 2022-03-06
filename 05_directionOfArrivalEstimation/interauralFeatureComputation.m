function [ipdCells, ivsCells, ildCells, AlgorithmStates, ipdDisambiguatedLogicalCells, itdCells, itdDisambiguatedLogicalCells] = ...
  interauralFeatureComputation(subbandSignalArray, AlgorithmParameters, ...
  AlgorithmStates)

    samplingRateHz = AlgorithmParameters.Gammatone.samplingRateHz;
    centerFreqsHz = AlgorithmStates.L.GammatoneStates.analyzer.center_frequencies_hz;
    cycleDurationSeconds = 1./centerFreqsHz;
    tau = AlgorithmParameters.nCyclesTau.*cycleDurationSeconds;
    % a = exp(-1./(samplingRateHz.*tau));
    
    nBands = AlgorithmParameters.Gammatone.nBands;

    ipdCells = cell(1, nBands);
    ipdDisambiguatedLogicalCells = cell(1, nBands); % for evaluation only
    ivsCells = cell(1, nBands);
    ildCells = cell(1, nBands);
    itdCells = cell(1, nBands);
    itdDisambiguatedLogicalCells = cell(1, nBands); % for evaluation only
    
    for iBand = 1:nBands
        subbandSignal.L = subbandSignalArray.L(:,iBand);
        subbandSignal.R = subbandSignalArray.R(:,iBand);
        States.L = AlgorithmStates.L.ProcessingStates{iBand};
        States.R = AlgorithmStates.R.ProcessingStates{iBand};
        States.Binaural = AlgorithmStates.Binaural.ProcessingStates{iBand};
        tauS = tau(iBand);
        
        %% ILD
        [subbandSignalLp, States] = ...
            ildBinauralLowPassFilter(States, subbandSignal, samplingRateHz);
        
        % % interaural level difference, eq. 5 in Dietz (2011)
        ildDb = 20*log10(abs(subbandSignalLp.L) ./ abs(subbandSignalLp.R));
        % % max(sig,1e-4) avoids division by zero
        %     20.*log10(max(subbandSignalLp.L,1e-4) ./ max(subbandSignalLp.L,1e-4));
        
        %% ITF
        itf = subbandSignal.L .* conj(subbandSignal.R);
        
        %% IVS
        [ivsMask, States.Binaural] = calcIvs(itf, tauS, samplingRateHz, ...
            AlgorithmParameters.ivsThreshold, States.Binaural);
%         ivsMask = true(size(ivsMask));

        %% IPD
        ipdRad = angle(itf);
        % test for equality: !!!
        % ipd_lp = angle(lowpass(outp.itf,a)); % Dietz2011 lowpass function
        [ipdRad, States.Binaural.ipdLP] = ...
            firstOrderLowPass(ipdRad, States.Binaural.ipdLP, AlgorithmParameters, tauS);
        
        ipdRadDisambiguated = disambiguateIpd(ipdRad, ildDb);
        ipdDisambiguatedLogical = ipdRad ~= ipdRadDisambiguated;

        %% ITD

        % instantaneous frequency (f_inst), eq. 4 in Dietz (2011)
        [instFreqLeft, States.L]  = calcInstFreq(subbandSignal.L,...
            samplingRateHz, States.L);
        [instFreqRight, States.R] = calcInstFreq(subbandSignal.R,...
            samplingRateHz, States.R);
        instFreq = max(eps,0.5*(instFreqLeft + instFreqRight)); % to avoid division by zero

        % interaural time difference (ITD), based on instantaneous frequencies
        itdSec = 1/(2*pi)*ipdRad./instFreq;

        itdSecDisambiguated = dietz2011_unwrapitd(itdSec,ildDb,instFreq);
        itdDisambiguatedLogical = itdSec ~= itdSecDisambiguated;
        %% write states back into global structs
        AlgorithmStates.L.ProcessingStates{iBand} = States.L;
        AlgorithmStates.R.ProcessingStates{iBand} = States.R;
        AlgorithmStates.Binaural.ProcessingStates{iBand} = States.Binaural;
        
        %% write computed values into structs
        ipdCells{iBand} = ipdRadDisambiguated;
        ipdDisambiguatedLogicalCells{iBand} = ipdDisambiguatedLogical;
        ivsCells{iBand} = ivsMask;
        ildCells{iBand} = ildDb;
        itdCells{iBand} = itdSecDisambiguated;
        itdDisambiguatedLogicalCells{iBand} = itdDisambiguatedLogical;
    end
end

function [f_inst, States] = calcInstFreq(sig, fs, States)
  % function f_inst = calc_f_inst(sig,fs);
  %
  % Calculates instantaneous frequency from a complex (analytical) signal
  % using first order differences
  %
  % input parameters:
  %   sig  : complex (analytical) input signal
  %   fs   : sampling frequency of sig
  %
  % output values:
  %   f_inst:   vector of estimated inst. frequency values
  %
  % copyright: Universitaet Oldenburg
  % author   : volker hohmann
  % date     : 12/2004
  %
  % adapted by Jeffrey Thomsen for sample or block-based processing 03/2022

  % handle signal samples for relay processing
  sig = [States.instFreqPreviousValue; sig];
  States.instFreqPreviousValue = sig(end);
  
  % original computation
  sig = sig./(abs(sig)+eps);

  f_inst = [sig(2:end).*conj(sig(1:end-1))];
  f_inst = angle(f_inst)/2/pi*fs;
end

function itd = dietz2011_unwrapitd(itd,ild,f_inst,tr)
    %% ===== Checking of input parameters ===================================
    nargmin = 3;
    nargmax = 4;
    narginchk(nargmin,nargmax);
    if nargin==3
        tr = 2.5;
    end
    %% ===== Calculation ====================================================
    itd = itd + ...
        round( ... % this will be -1,0,1
            0.4*sign(round(ild/2 / (abs(tr)+1e-9))) - 0.4*sign(itd) ) ...
        ./ f_inst;
end