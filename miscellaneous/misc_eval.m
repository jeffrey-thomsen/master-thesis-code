% compute glimpse and ivs mask ratios

sumtarget = 0;for iBand= 1:30;
sumtarget = sumtarget+numel(SimulationData{6,4}.Data.targetSampleIndices.R{iBand});
end;
suminterf = 0;for iBand= 1:30;
suminterf = suminterf+numel(SimulationData{6,4}.Data.interfSampleIndices.R{iBand});
end;
glimpseRatio_J = (suminterf+sumtarget)/(30*36563)

sumivsmask = 0;for iBand = 1:30;
sumivsmask = sumivsmask+nnz(SimulationData{26,4}.Data.ivsMaskCells{iBand});
end;
sumivsmask/(30*36563)

% instationarity

nnz(diff(maskTarget_J.L',1,2))/(30*36563)
nnz(diff(maskTarget_J.R',1,2))/(30*36563)
nnz(diff(maskInterf_J.L',1,2))/(30*36563)
nnz(diff(maskInterf_J.R',1,2))/(30*36563)

%% IBM instationarity
iSp = 6; jAng = 4;
nnz(diff(ibmTarget{iSp,jAng}.L',1,2))/(30*36563)
nnz(diff(ibmTarget{iSp,jAng}.R',1,2))/(30*36563)
nnz(diff(ibmTarget{iSp,jAng}.L',1,2))/(30*36563)
nnz(diff(ibmTarget{iSp,jAng}.R',1,2))/(30*36563)

%% Show how glimpses vary with angle configurations
iSp = 6;
for jAng = 1:6
plotIbmGlimpses(maskTarget{iSp,jAng}, maskInterf{iSp,jAng}, ...
    ibmTarget{iSp,jAng}, ibmInterf{iSp,jAng}, inputTargetSignal{iSp,jAng}, ...
    inputInterfSignal{iSp,jAng}, SimulationData{iSp,jAng}, timeVec, ...
    AlgorithmParameters.Gammatone.samplingRateHz);
end