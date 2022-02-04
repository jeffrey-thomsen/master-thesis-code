function [sigma, delta] = calcSigmaDelta(subbandSignals, p0Values)
    sigma.L = subbandSignals.L; %+ subbandSignals.L(-p0Values);
    sigma.R = subbandSignals.R; %+ subbandSignals.R(-p0Values);
    delta.L = subbandSignals.L; %+ subbandSignals.L(-p0Values);
    delta.R = subbandSignals.R; %+ subbandSignals.R(-p0Values);
end