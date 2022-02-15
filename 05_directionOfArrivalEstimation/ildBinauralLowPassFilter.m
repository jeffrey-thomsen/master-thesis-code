function [subbandSignalLp, States] = ...
    ildBinauralLowPassFilter(States, subbandSignal, samplingRateHz)
    
% low pass filter with a fixed cutoff frequency for every frequency channel
    [b,a] = butter(2,30/(samplingRateHz/2),'low');
    
    [subbandSignalLp.L, States.L.ildLP] = ...
        filter(b, a, subbandSignal.L, States.L.ildLP);
   
    [subbandSignalLp.R, States.R.ildLP] = ...
        filter(b, a, subbandSignal.R, States.R.ildLP);

%     when is this identical to my simple lowpass?? tau=1/fc or
%     1/(2*pi*fc)???
%     [subbandSignalLp.L, States.L.ildLP] = ...
%         firstOrderLowPass(subbandSignal.L, States.L.ildLP, AlgorithmParameters, 1/30);
%     [subbandSignalLp.R, States.R.ildLP] = ...
%         firstOrderLowPass(subbandSignal.R, States.R.ildLP, AlgorithmParameters, 1/30);
end