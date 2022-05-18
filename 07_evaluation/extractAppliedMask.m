% Construct logical matrices from the time-frequency target and interferer
% glimpses stored in the SimulationData struct for comparing them to IBM
function [maskTarget, maskInterf] = extractAppliedMask(SimulationData, nBands)

    signalIndices = SimulationData.signalIndices;
    for iBand = 1:nBands
            maskTarget.L(:,iBand) = ismember(signalIndices, ...
                SimulationData.Data.targetSampleIndices.L{iBand});
            maskTarget.R(:,iBand) = ismember(signalIndices, ...
                SimulationData.Data.targetSampleIndices.R{iBand});
            maskInterf.L(:,iBand) = ismember(signalIndices, ...
                SimulationData.Data.interfSampleIndices.L{iBand});
            maskInterf.R(:,iBand) = ismember(signalIndices, ...
                SimulationData.Data.interfSampleIndices.R{iBand});
    end
end