function [] = plotIbmGlimpses(maskTarget, maskInterf, ibmTarget, ibmInterf, ...
  targetSignal, interfSignal, SimulationData, timeVec, samplingRateHz)    

    p0SearchRangeSecondsVector = ...
        SimulationData.Data.p0SearchRangeSamplesVector ./ samplingRateHz;
    p0SearchRangeFreqHzVector = 1./p0SearchRangeSecondsVector;


    figure;
    plotMonoSpectrogram(targetSignal, samplingRateHz);
    ylim([0, 0.5]);
%     caxis([-80 -10]);

    hold on;
    for iBand = 1:30
        tmpP0IndexVector = SimulationData.Data.p0DetectedIndexVectors.L{iBand};

        tmpIndices = (maskTarget.L(:,iBand) & ibmTarget.L(:,iBand));
        scatter(timeVec(tmpIndices), ...
                1e-3.*p0SearchRangeFreqHzVector(tmpP0IndexVector(tmpIndices)), ...
                15, 'green','filled');

        tmpIndices = (maskTarget.L(:,iBand) & ~ibmTarget.L(:,iBand));
        scatter(timeVec(tmpIndices), ...
                1e-3.*p0SearchRangeFreqHzVector(tmpP0IndexVector(tmpIndices)), ...
                15, 'red','filled');

        
        tmpP0IndexVector = SimulationData.Data.p0DetectedIndexVectors.R{iBand};

        tmpIndices = (maskTarget.R(:,iBand) & ibmTarget.R(:,iBand));
        scatter(timeVec(tmpIndices), ...
                1e-3.*p0SearchRangeFreqHzVector(tmpP0IndexVector(tmpIndices)), ...
                15, 'green','filled');
        
        tmpIndices = (maskTarget.R(:,iBand) & ~ibmTarget.R(:,iBand));
        scatter(timeVec(tmpIndices), ...
                1e-3.*p0SearchRangeFreqHzVector(tmpP0IndexVector(tmpIndices)), ...
                15, 'red','filled');
    end
    legend('correct target glimpses', 'incorrect target glimpses')
    title('Target mask')

    figure;

    plotMonoSpectrogram(interfSignal, samplingRateHz);
    ylim([0, 0.5]);
%     caxis([-80 -10]);

    hold on;
    for iBand = 1:30
        tmpP0IndexVector = SimulationData.Data.p0DetectedIndexVectors.L{iBand};

        tmpIndices = (maskInterf.L(:,iBand) & ibmInterf.L(:,iBand));
        scatter(timeVec(tmpIndices), ...
                1e-3.*p0SearchRangeFreqHzVector(tmpP0IndexVector(tmpIndices)), ...
                15, 'green','filled');

        tmpIndices = (maskInterf.L(:,iBand) & ~ibmInterf.L(:,iBand));
        scatter(timeVec(tmpIndices), ...
                1e-3.*p0SearchRangeFreqHzVector(tmpP0IndexVector(tmpIndices)), ...
                15, 'red','filled');

        
        tmpP0IndexVector = SimulationData.Data.p0DetectedIndexVectors.R{iBand};

        tmpIndices = (maskInterf.R(:,iBand) & ibmInterf.R(:,iBand));
        scatter(timeVec(tmpIndices), ...
                1e-3.*p0SearchRangeFreqHzVector(tmpP0IndexVector(tmpIndices)), ...
                15, 'green','filled');
        
        tmpIndices = (maskInterf.R(:,iBand) & ~ibmInterf.R(:,iBand));
        scatter(timeVec(tmpIndices), ...
                1e-3.*p0SearchRangeFreqHzVector(tmpP0IndexVector(tmpIndices)), ...
                15, 'red','filled');
    end
    legend('correct interferer glimpses', 'incorrect interferer glimpses')
    title('Interferer mask')

end