%% visualize target and interferer glimpses
function [] = plotPeriodGlimpses(SimulationData, timeVec, samplingRateHz)
for iBand = [3,10,25]
    figure;
    title('Gammatone band no.',num2str(iBand))
    
    subplot(1,2,1)
    title('Processed periodic samples - Left channel')
    xlabel('time(s)')
    ylabel('detected frequency (Hz)')
    hold on;
    scatter(timeVec(find(SimulationData.p0DetectedIndexVectors.L{iBand})), ...
        samplingRateHz ./ SimulationData.p0SearchRangeSamplesVector(...
        nonzeros(SimulationData.p0DetectedIndexVectors.L{iBand})),...
        15, 'green', 'filled');
    legend('glimpses')
    ylim([100, 350])
    
    subplot(1,2,2)
    title('Processed periodic samples - Right channel')
    xlabel('time (s)')
    ylabel('detected frequency (Hz)')
    hold on;
    scatter(timeVec(find(SimulationData.p0DetectedIndexVectors.R{iBand})), ...
        samplingRateHz ./ SimulationData.p0SearchRangeSamplesVector(...
        nonzeros(SimulationData.p0DetectedIndexVectors.R{iBand})),...
        15, 'green', 'filled');
    legend('glimpses')
    ylim([100, 350])
end