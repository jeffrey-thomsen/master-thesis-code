figure;
iPlot = 0;
for iAn = [1:4:20, 2:4:20, 3:4:20, 4:4:20]
    iPlot = iPlot + 1;
    azDeg = extractAzDeg(SimulationData{1,iAn}.Data.ivsMaskCells, SimulationData{1,iAn}.Data.azimuthDegCells, AlgorithmParameters, 30);
    for iSp=2:45
        azDeg = [azDeg;extractAzDeg(SimulationData{iSp,iAn}.Data.ivsMaskCells, SimulationData{iSp,iAn}.Data.azimuthDegCells, AlgorithmParameters, 30)];
    end
    subplot(4,5,iPlot)
    doaHistogram(azDeg)
    title([num2str(anglePermutations(iAn,1)),' ',num2str(anglePermutations(iAn,2))])
end