function [resynthesizedSignal, synthesizer] = ...
  subbandResynthesis(subbandSignals, synthesizer)
    
    [resynthesizedSignal, synthesizer] = ...
        hohmann2002_process(synthesizer, subbandSignals);
    
end