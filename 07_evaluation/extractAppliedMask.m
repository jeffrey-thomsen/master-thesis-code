% Construct logical matrices from the time-frequency target and interferer
% glimpses stored in the SimulationData struct for comparing them to IBM
function [maskTarget, maskInterf] = extractAppliedMask(SimulationData, nBands)

%     targetSampleIndices.L = cell(size(SimulationData));
%     targetSampleIndices.R = cell(size(SimulationData));
%     interfSampleIndices.L = cell(size(SimulationData));
%     interfSampleIndices.R = cell(size(SimulationData));
%     maskTarget = cell(size(SimulationData));
%     maskInterf = cell(size(SimulationData));

%     for iSignal = 1:numel(SimulationData)
        signalIndices = SimulationData.signalIndices;
        targetSampleIndices.L = cell(1, nBands);
        targetSampleIndices.R = cell(1, nBands);
        interfSampleIndices.L = cell(1, nBands);
        interfSampleIndices.R = cell(1, nBands);
        for iBand = 1:nBands
            for iBlock = 1:numel(SimulationData.Data)
                tmpTargetIndicesReBlock.L = ...
                    SimulationData.Data(iBlock).targetSampleIndices.L{iBand};
                tmpTargetIndicesReSignal.L = signalIndices(iBlock, tmpTargetIndicesReBlock.L);
                targetSampleIndices.L{iBand} = ...
                    [targetSampleIndices.L{iBand}, ...
                    tmpTargetIndicesReSignal.L];
                               
                tmpTargetIndicesReBlock.R = ...
                    SimulationData.Data(iBlock).targetSampleIndices.R{iBand};
                tmpTargetIndicesReSignal.R = signalIndices(iBlock, tmpTargetIndicesReBlock.R);
                targetSampleIndices.R{iBand} = ...
                    [targetSampleIndices.R{iBand}, ...
                    tmpTargetIndicesReSignal.R];
                                
                tmpInterfIndicesReBlock.L = ...
                    SimulationData.Data(iBlock).interfSampleIndices.L{iBand};
                tmpInterfIndicesReSignal.L = signalIndices(iBlock, tmpInterfIndicesReBlock.L);
                interfSampleIndices.L{iBand} = ...
                    [interfSampleIndices.L{iBand}, ...
                    tmpInterfIndicesReSignal.L];
                                
                tmpInterfIndicesReBlock.R = ...
                    SimulationData.Data(iBlock).interfSampleIndices.R{iBand};
                tmpInterfIndicesReSignal.R = signalIndices(iBlock, tmpInterfIndicesReBlock.R);
                interfSampleIndices.R{iBand} = ...
                    [interfSampleIndices.R{iBand}, ...
                    tmpInterfIndicesReSignal.R];
                
            end
        end
%     end

    nSamples = max(max(signalIndices));
%     for iSignal = 1:numel(SimulationData)
        for iBand = 1:nBands
            maskTarget.L(:,iBand) = ismember(1:nSamples, ...
                targetSampleIndices.L{iBand});
            maskTarget.R(:,iBand) = ismember(1:nSamples, ...
                targetSampleIndices.R{iBand});
            maskInterf.L(:,iBand) = ismember(1:nSamples, ...
                interfSampleIndices.L{iBand});
            maskInterf.R(:,iBand) = ismember(1:nSamples, ...
                interfSampleIndices.R{iBand});
        end
%     end
end