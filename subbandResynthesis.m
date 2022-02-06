% Runs the Hohmann2002 gammatone synthesis filterbank
% subbandSignals - complex-valued set of subband signals to be
% resynthesized
% synthesizer - struct containing the filterbank parameters and states
% resynthesizedSignal - Nx1 real-valued signal
function [resynthesizedSignal, synthesizer] = ...
  subbandResynthesis(subbandSignals, synthesizer)
    
    [resynthesizedSignal, synthesizer] = ...
        hohmann2002_process(synthesizer, subbandSignals);
    
end