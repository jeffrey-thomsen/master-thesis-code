% Interaural feature computaiton of the Thomsen2022 speech enhancement
% algorithm. Uses a set of binaural signals that has been decomposed by a
% gammatone filterbank and returns IPD, ILD and ITD values and an IVS
% coherence mask for each subband sample
%
% N - length of signal
% M - number of subbands
%
% Input:
% subbandSignalArray - struct containing NxM arrays of complex-valued
% subband samples from a gammatone filterbank, for left and right channel
% AlgorithmParameters - struct containing all simulation parameters
% AlgorithmStates - struct containing filter states and FIFO arrays that
% need to be handed over to the next signal block/sample to be processed
%
% Output:
% ipdCells, ildCells, itdCells - cells of size M containing IPD values in
% rad, ILD values in dB and ITD values in seconds for each subband sample
% ivsCells - cell of size M containing the coherence masks for each subband
% signal in the form of a logical vector
% ipdDisambiguatedLogicalCells, itdDisambiguatedLogicalCells - cells of
% size M containing logical arrays determining whether an IPD or ITD value
% was disambiguated by the ILD criterion (for algorithm evaluation purposes
% only)
function [ipdCells, ivsCells, ildCells, AlgorithmStates, ipdDisambiguatedLogicalCells, itdCells, itdDisambiguatedLogicalCells] = ...
  interauralFeatureComputation(subbandSignalArray, AlgorithmParameters, ...
  AlgorithmStates)

    samplingRateHz = AlgorithmParameters.Gammatone.samplingRateHz;

    % calculate time constants for the LP filters used in each subband
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

function [instFrequency, States] = ...
  calcInstFreq(signal, samplingRateHz, States)
    % function f_inst = calcInstFreq(sig,fs);
    %
    % Calculates instantaneous frequency from a complex (analytical) signal
    % using first order differences
    %
    % input parameters:
    %   signal: complex (analytical) input signal
    %   samplingRateHz: sampling frequency of signal
    %   States: struct containing the final signal value from the
    %   previously processed block
    %
    % output values:
    %   instFrequency: vector of estimated instantaneous frequency values
    %   States: see above
    %
    % copyright: Universitaet Oldenburg
    % author   : volker hohmann
    % date     : 12/2004
    %
    % adapted by Jeffrey Thomsen for sample or block-based processing 03/2022
    
    % handle signal samples for relay processing
    signal = [States.instFreqPreviousValue; signal];
    States.instFreqPreviousValue = signal(end);
    
    % original computation
    signal = signal./(abs(signal)+eps);
    
    instFrequency = signal(2:end) .* conj(signal(1:end-1));
    instFrequency = angle(instFrequency)/2/pi*samplingRateHz;
end

function itd = dietz2011_unwrapitd(itd, ild, instFrequency, thresholdDb)
    % author: Mathias Dietz
    % set the ITD to the sign of the ILD if it is above a certain threshold
    %% Checking of input parameters
    nargmin = 3;
    nargmax = 4;
    narginchk(nargmin,nargmax);
    if nargin==3
        thresholdDb = 2.5;
    end
    %% Calculation
    itd = itd + ...
        round( ... % this will be -1,0,1
            0.4*sign(round(ild/2 / (abs(thresholdDb)+1e-9))) - 0.4*sign(itd) ) ...
        ./ instFrequency;
end