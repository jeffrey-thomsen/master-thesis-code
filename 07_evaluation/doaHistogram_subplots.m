% script for generating the DOA histogram subplots per DOA angle combination
figure;
iPlot = 0;

if nAnglePerms == 20
    nAnPlots = [1:4:20, 2:4:20, 3:4:20, 4:4:20];
elseif nAnglePerms == 6
    nAnPlots = [1,3,5,2,4,6];
end

for iAn = nAnPlots
    iPlot = iPlot + 1;
    if nSpeakerCombos == 22
        iStart = 24;
    else
        iStart = 1;
    end
    azDeg = extractAzDeg(SimulationData{iStart,iAn}.Data.ivsMaskCells, SimulationData{iStart,iAn}.Data.azimuthDegCells, AlgorithmParameters, 30);
    for iSp=iStart+1:nSpeakerCombos
        azDeg = [azDeg;extractAzDeg(SimulationData{iSp,iAn}.Data.ivsMaskCells, SimulationData{iSp,iAn}.Data.azimuthDegCells, AlgorithmParameters, 30)];
    end

    if nAnglePerms == 20
        subplot(4,5,iPlot)
    elseif nAnglePerms == 6
        subplot(2,3,iPlot)
    end

    doaHistogram(azDeg)
    title(['DOA ID: ',num2str(iAn),' - $\varphi_{0,t}$: ',num2str(anglePermutations(iAn,1)),'$^{\circ}$, $\varphi_{0,i}$: ',num2str(anglePermutations(iAn,2)),'$^{\circ}$'],'Interpreter','latex')
    if iAn<5
        ylabel('Probability of occurrence')
    end
    if mod(iAn,4) == 0
        xlabel('$\hat{\varphi} (^{\circ})$','Interpreter','latex')
    end

    xlim([-95 95])
    ylim([0 0.2])
    clear azDeg
end