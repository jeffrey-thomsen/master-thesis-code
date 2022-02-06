% Takes a binaural pair of subband signals, corresponding gammatone
% synthesis filterbank paramters, filter states and hands the individual
% left and right channels to subbandResynthesis.m
% subbandSignals - struct containing a set of complex-valued subband 
% signals for the left and right channel
% AlgorithmParameters - struct containing a struct with filter parameters
% and filter states for the gammatone synthesis filterbank for each the 
% left and right channel signals
% resynthesizedSignal - Nx2 matrix containing the real-valued binaural signal
function [resynthesizedSignal, AlgorithmParameters] = ...
  subbandResynthesisBinaural(subbandSignals, AlgorithmParameters)
    
    % read the filter parameters and states for left and right channel
    % synthesis filterbanks
    synthesizer.L = AlgorithmParameters.L.FilterStates.Gammatone.synthesizer;
    synthesizer.R = AlgorithmParameters.R.FilterStates.Gammatone.synthesizer;
    
    % resynthesize left and right channel binaural signals
    resynthesizedSignal = zeros(length(subbandSignals.L),2);
    
    [resynthesizedSignal(:,1), synthesizer.L] = ...
        subbandResynthesis(subbandSignals.L, synthesizer.L);
    
    [resynthesizedSignal(:,2), synthesizer.R] = ...
        subbandResynthesis(subbandSignals.R, synthesizer.R);
    
    % update filter parameters
    AlgorithmParameters.L.FilterStates.Gammatone.synthesizer = synthesizer.L;
    AlgorithmParameters.R.FilterStates.Gammatone.synthesizer = synthesizer.R;
    
end