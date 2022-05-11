% Construct logical matrices from the time-frequency target and interferer
% glimpses stored in the SimulationData struct for comparing them to IBM
function [maskTarget, maskInterf] = extractAppliedMask(SimulationData, nBands)

    targetSampleIndices.L = cell(size(SimulationData));
    targetSampleIndices.R = cell(size(SimulationData));
    interfSampleIndices.L = cell(size(SimulationData));
    interfSampleIndices.R = cell(size(SimulationData));
    maskTarget = cell(size(SimulationData));
    maskInterf = cell(size(SimulationData));

    for iSignal = 1:numel(SimulationData)
        signalIndices = SimulationData{iSignal}.signalIndices;
        targetSampleIndices.L{iSignal} = cell(1, nBands);
        targetSampleIndices.R{iSignal} = cell(1, nBands);
        interfSampleIndices.L{iSignal} = cell(1, nBands);
        interfSampleIndices.R{iSignal} = cell(1, nBands);
        for iBand = 1:nBands
            for iBlock = 1:numel(SimulationData{iSignal}.Data)
                tmpTargetIndicesReBlock.L = ...
                    SimulationData{iSignal}.Data(iBlock).targetSampleIndices.L{iBand};
                tmpTargetIndicesReSignal.L = signalIndices(iBlock, tmpTargetIndicesReBlock.L);
                targetSampleIndices.L{iSignal}{iBand} = ...
                    [targetSampleIndices.L{iSignal}{iBand}, ...
                    tmpTargetIndicesReSignal.L];
                               
                tmpTargetIndicesReBlock.R = ...
                    SimulationData{iSignal}.Data(iBlock).targetSampleIndices.R{iBand};
                tmpTargetIndicesReSignal.R = signalIndices(iBlock, tmpTargetIndicesReBlock.R);
                targetSampleIndices.R{iSignal}{iBand} = ...
                    [targetSampleIndices.R{iSignal}{iBand}, ...
                    tmpTargetIndicesReSignal.R];
                                
                tmpInterfIndicesReBlock.L = ...
                    SimulationData{iSignal}.Data(iBlock).interfSampleIndices.L{iBand};
                tmpInterfIndicesReSignal.L = signalIndices(iBlock, tmpInterfIndicesReBlock.L);
                interfSampleIndices.L{iSignal}{iBand} = ...
                    [interfSampleIndices.L{iSignal}{iBand}, ...
                    tmpInterfIndicesReSignal.L];
                                
                tmpInterfIndicesReBlock.R = ...
                    SimulationData{iSignal}.Data(iBlock).interfSampleIndices.R{iBand};
                tmpInterfIndicesReSignal.R = signalIndices(iBlock, tmpInterfIndicesReBlock.R);
                interfSampleIndices.R{iSignal}{iBand} = ...
                    [interfSampleIndices.R{iSignal}{iBand}, ...
                    tmpInterfIndicesReSignal.R];
                
            end
        end
    end

    nSamples = max(max(signalIndices));
    for iSignal = 1:numel(SimulationData)
        for iBand = 1:nBands
            maskTarget{iSignal}.L(:,iBand) = ismember(1:nSamples, ...
                targetSampleIndices.L{iSignal}{iBand});
            maskTarget{iSignal}.R(:,iBand) = ismember(1:nSamples, ...
                targetSampleIndices.R{iSignal}{iBand});
            maskInterf{iSignal}.L(:,iBand) = ismember(1:nSamples, ...
                interfSampleIndices.L{iSignal}{iBand});
            maskInterf{iSignal}.R(:,iBand) = ismember(1:nSamples, ...
                interfSampleIndices.R{iSignal}{iBand});
        end
    end
end