% Find out which of the 5 HRTF measurements is "most median", i.e. which
% measurement is the median amongst the 5 for the majority of all HRIR
% samples

%% load impulse responses
sofaStructs = cell(1,5);
sofaStructs{1} = SOFAload('HRIR_KEMAR_DV0001_1.sofa',[5 2],'R');
sofaStructs{2} = SOFAload('HRIR_KEMAR_DV0001_2.sofa',[5 2],'R');
sofaStructs{3} = SOFAload('HRIR_KEMAR_DV0001_3.sofa',[5 2],'R');
sofaStructs{4} = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
sofaStructs{5} = SOFAload('HRIR_KEMAR_DV0001_5.sofa',[5 2],'R');

hrir = cell(1,5);
for i=1:5
    hrir{i} = sofaStructs{i}.Data.IR;
end

%% iterate over all incidence angles and both ear channels
medianCounts = zeros(5,1);
for jAngle = 1:87
    for kChannel = 1:2

        % load the impulse response of each measurement
        tempHRIR = zeros(356,5);
        for iMeasurement = 1:5
            tempHRIR(:,iMeasurement) = ...
                squeeze(hrir{iMeasurement}(jAngle,kChannel,:));
        end
        % compute the median value for each sample of the impulse responses
        tempMedian = median(tempHRIR,2);
        % for each measurement, count the number of samples that are the 
        % median value
        for iMeasurement = 1:5
            medianCounts(iMeasurement) = medianCounts(iMeasurement) + ...
                numel(nonzeros(tempHRIR(:,iMeasurement)==tempMedian));
        end

    end
end

%% show result
[maxVal, maxIndex] = max(medianCounts);
disp(['The best HRIR to use is HRIR_KEMAR_DV0001_', num2str(maxIndex),'.sofa'])
disp(['It contains ', num2str(100*((maxVal/mean(medianCounts))-1)), '% more median counts than the mean of median counts among the measurements.'])