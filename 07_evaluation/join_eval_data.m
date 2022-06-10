tic

load('2022-05-31_10-36_simulation_data.mat');

fprintf('\n')
for iSp = 24:45
    fprintf('%i,',iSp)
    for iAn = 1:20
        SimulationData{iSp,iAn} = rmfield(SimulationData{iSp,iAn}, 'signalIndices');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'p0SearchRangeSamplesVector');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'p0DetectedIndexVectors');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'targetSampleIndices');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'interfSampleIndices');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'snrDesired');
%         SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'azimuthDegCells');
    end
end

SimulationData2 = SimulationData;
clear SimulationData

load('2022-05-30_16-48_simulation_data.mat');

fprintf('\n')
for iSp = 1:23
    fprintf('%i,',iSp)
    for iAn = 1:20
        SimulationData{iSp,iAn} = rmfield(SimulationData{iSp,iAn}, 'signalIndices');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'p0SearchRangeSamplesVector');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'p0DetectedIndexVectors');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'targetSampleIndices');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'interfSampleIndices');
        SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'snrDesired');
%         SimulationData{iSp,iAn}.Data = rmfield(SimulationData{iSp,iAn}.Data, 'azimuthDegCells');
    end
end

fprintf('\n')
for iSp = 24:45
    fprintf('%i,',iSp)
    for iAn = 1:20
        SimulationData{iSp,iAn} = SimulationData2{iSp,iAn};
    end
end

% SimulationData2 = SimulationData;
% clear SimulationData
% 
% load('2022-05-29_23-03_azDeg_simulation_data.mat');
% 
% fprintf('\n')
% for iSp = 1:45
%     fprintf('%i,',iSp)
%     for iAn = 1:20
%         SimulationData{iSp,iAn}.Data.ivsMaskCells = SimulationData2{iSp,iAn}.Data.ivsMaskCells;
%     end
% end

save('2022-06-02_13-05-03_DOA_simulation_data.mat', ...
    'AlgorithmParameters', 'TestSignalParameters', ...
    'anglePermutations', 'speakerCombinations', ...
    'SimulationData', '-v7.3');

toc
