function [resynthesizedSignal, AlgorithmParameters] = ...
  subbandResynthesisBinaural(subbandSignals, AlgorithmParameters)

    synthesizer.L = AlgorithmParameters.L.FilterStates.Gammatone.synthesizer;
    synthesizer.R = AlgorithmParameters.R.FilterStates.Gammatone.synthesizer;
    
    resynthesizedSignal = zeros(length(subbandSignals.L),2);
    [resynthesizedSignal(:,1), synthesizer.L] = ...
        subbandResynthesis(subbandSignals.L, synthesizer.L);
    [resynthesizedSignal(:,2), synthesizer.R] = ...
        subbandResynthesis(subbandSignals.R, synthesizer.R);

    AlgorithmParameters.L.FilterStates.Gammatone.synthesizer = synthesizer.L;
    AlgorithmParameters.R.FilterStates.Gammatone.synthesizer = synthesizer.R;
    
end