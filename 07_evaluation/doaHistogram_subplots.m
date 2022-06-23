% script for generating the DOA histogram subplots per DOA angle combination

iPlot = 0;

if nAnglePerms == 20
    figure('Position',[1495 40 1065 1300]);
    nAnPlots = 1:20;%[1:4:20, 2:4:20, 3:4:20, 4:4:20];
    subplot = @(m,n,p) subtightplot(m,n,p,[0.05 0.03], [0.06 0.03], [0.06 0.02]);
elseif nAnglePerms == 6
    figure('Position',[1500 800 1065 540]);
    nAnPlots = [1,3,5,2,4,6];
    subplot = @(m,n,p) subtightplot(m,n,p,[0.11 0.05], [0.1 0.05], [0.06 0.02]);
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
        subplot(5,4,iPlot)
    elseif nAnglePerms == 6
        subplot(2,3,iPlot)
    end

    doaHistogram(azDeg)
    title(['DOA ID:',num2str(iAn),' - $\varphi_{0,t}$:',num2str(anglePermutations(iAn,1)),'$^{\circ}$, $\varphi_{0,i}$:',num2str(anglePermutations(iAn,2)),'$^{\circ}$'],'Interpreter','latex')
    
    if nAnglePerms == 20
        if ismember(iAn,1:4:20)%iAn<5
            ylabel('Probability of occurrence')
        else
            set(gca,'YTickLabel',[]);
        end
        if iAn>16%mod(iAn,4) == 0
            xlabel('$\hat{\varphi} (^{\circ})$','Interpreter','latex')
        end
    elseif nAnglePerms == 6
        if ismember(iAn,[1,2])
            ylabel('Probability of occurrence')
        end
        if ismember(iAn,[2,4,6])
            xlabel('$\hat{\varphi} (^{\circ})$','Interpreter','latex')
        end
    end

    xlim([-95 95])
    ylim([0 0.2])
    clear azDeg
end