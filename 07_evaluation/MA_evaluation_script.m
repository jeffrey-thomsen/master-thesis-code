% script for generating MA plots
clear

%% Fetch DOA Data

%% DOA histograms general
% tic
% 
% % dateString = '2022-05-27_02-55_DOA'; % CFR0, IVS0, full DOA
% % dateString = '2022-05-27_21-16'; % CFR0, full DOA
% 
% % dateString = '2022-05-28_16-04'; nSpeakerCombos = 23; % full DOA speakers 1:23
% % dateString = '2022-05-29_09-40'; nSpeakerCombos = 22; % full DOA speakers 24:45
% % dateString = '2022-05-29_23-03_DOA'; % full DOA combined
% 
% % dateString = '2022-05-30_16-48'; nSpeakerCombos = 23; % simple DOA speakers 1:23
% % dateString = '2022-05-31_10-36'; nSpeakerCombos = 22; % simple DOA speakers 24:45
% % dateString = '2022-06-02_13-05_DOA'; % simple DOA combined
% 
% % Control set
% % dateString = '2022-06-03_15-17'; % CFR0, IVS0, full DOA
% % dateString = '2022-06-03_20-08'; % CFR0, full DOA
% % dateString = '2022-06-05_02-48'; % Full DOA
% % dateString = '2022-06-05_11-33'; % Random p0
% % dateString = '2022-06-05_14-41'; % ChenP0Detection
% % dateString = '2022-06-05_22-59'; % AzimuthPooling
% % dateString = '2022-06-06_03-11'; % Enhancement
% % dateString = '2022-06-06_11-04'; % Cancellation
% 
% load([dateString,'_simulation_data.mat']);
% 
% nSpeakerCombos = size(speakerCombinations,1);
% nAnglePerms = size(anglePermutations,1);
% 
% doaHistogram_subplots;
% 
% toc

%% DOA histograms specific

tic
load('2022-05-27_02-55_DOA_simulation_data.mat');
nSpeakerCombos = size(speakerCombinations,1);
nAnglePerms = size(anglePermutations,1);
doaHistogram_subplots;
clear SimulationData
toc

tic
load('2022-05-29_23-03_DOA_simulation_data.mat');
nSpeakerCombos = size(speakerCombinations,1);
nAnglePerms = size(anglePermutations,1);
doaHistogram_subplots;
toc

tic
SimulationData = SimulationData(1:2:45,:);
SimulationData = SimulationData(:,[2,3,9,11,13,15]);
speakerCombinations = speakerCombinations(1:2:45,:);
anglePermutations = anglePermutations([2,3,9,11,13,15],:);
nSpeakerCombos = size(speakerCombinations,1);
nAnglePerms = size(anglePermutations,1);
doaHistogram_subplots;
clear SimulationData
toc

tic
load('2022-06-05_02-48_simulation_data.mat');
nSpeakerCombos = size(speakerCombinations,1);
nAnglePerms = size(anglePermutations,1);
doaHistogram_subplots;
clear SimulationData
toc

%% Fetch data for rest of plots

% tic
% 
% % dateString = '2022-05-27_02-55'; % CFR0, IVS0, full DOA
% % dateString = '2022-05-27_21-16'; % CFR0, full DOA
% 
% % dateString = '2022-05-29_23-03'; % full DOA
% 
% % dateString = '2022-06-02_13-05'; % simple DOA
% 
% % control set
% dateString = '2022-06-03_15-17'; % CFR0, IVS0, full DOA
% % dateString = '2022-06-03_20-08'; % CFR0, full DOA
% % dateString = '2022-06-05_02-48'; % Full DOA
% % dateString = '2022-06-05_11-33'; % Random p0
% % dateString = '2022-06-05_14-41'; % ChenP0Detection
% % dateString = '2022-06-05_22-59'; % AzimuthPooling
% % dateString = '2022-06-06_03-11'; % Enhancement
% % dateString = '2022-06-06_11-04'; % Cancellation
% 
% load([dateString,'_evaluation_data.mat']);
% % load([dateString,'_signal_data.mat']);
% 
% clear ibmTarget
% clear ibmInterf
% clear maskAppliedTarget
% clear maskAppliedInterf
% 
% % Load and plot IBM+TIR data
% 
% [precisionTarget, precisionInterf, recallTarget, recallInterf, ...
%     F1Target, F1Interf] = extractIbmEvalVals(speakerCombinations, ...
%     anglePermutations, precision, recall, deltaSNR_H);
% 
% filename = [dateString,'_meta_eval_data.mat'];
% save(filename, ...
% 'deltaSNR_H', 'F1Target', 'F1Interf',...
% 'precisionTarget', 'precisionInterf',...
% 'recallTarget', 'recallInterf',...
% 'MetaData', '-v7.3');
% 
% toc
% 
% % Display IBM+TIR means+stds
% 
% displayMeanStd(precisionTarget,recallTarget,precisionInterf,recallInterf,F1Target,F1Interf,deltaSNR_H)



%% Load previous meta evaluation data for plotting

ThreshYes = load('2022-05-27_02-55_meta_eval_data.mat');
displayMeanStd(ThreshYes.precisionTarget, ThreshYes.recallTarget, ...
    ThreshYes.precisionInterf, ThreshYes.recallInterf, ...
    ThreshYes.F1Target, ThreshYes.F1Interf, ThreshYes.deltaSNR_H)
ThreshNo = load('2022-05-29_23-03_meta_eval_data.mat');
displayMeanStd(ThreshNo.precisionTarget, ThreshNo.recallTarget,...
    ThreshNo.precisionInterf, ThreshNo.recallInterf, ...
    ThreshNo.F1Target, ThreshNo.F1Interf, ThreshNo.deltaSNR_H)
%% Precision+recall Hisograms+Heatmaps
plotDOAHistHeat(ThreshYes.precisionTarget,ThreshYes.precisionInterf,ThreshYes.recallTarget,ThreshYes.recallInterf)
plotDOAHistHeat(ThreshNo.precisionTarget,ThreshNo.precisionInterf,ThreshNo.recallTarget,ThreshNo.recallInterf)

%% Comparison boxplots
figure;
hold on
% boxchart(-0.3*ones(900,1),ThreshYes.F1Target(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
% boxchart(0.3*ones(900,1),ThreshNo.F1Target(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(1.7*ones(900,1),ThreshYes.precisionTarget(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(2.3*ones(900,1),ThreshNo.precisionTarget(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(3.7*ones(900,1),ThreshYes.recallTarget(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(4.3*ones(900,1),ThreshNo.recallTarget(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
% boxchart(5.7*ones(900,1),ThreshYes.F1Interf(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
% boxchart(6.3*ones(900,1),ThreshNo.F1Interf(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(7.7*ones(900,1),ThreshYes.precisionInterf(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(8.3*ones(900,1),ThreshNo.precisionInterf(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(9.7*ones(900,1),ThreshYes.recallInterf(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(10.3*ones(900,1),ThreshNo.recallInterf(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
legend('with thresholding','without thresholding')
xticks([2, 4, 8, 10])
xticklabels({'Target precision', 'Target recall', 'Interferer precision', 'Interferer recall'})
ylabel('Rate (1)')

%
% figure;
% hold on
% yyaxis left
% ylabel('Rate (1)')
% boxchart(0*ones(900,1),ThreshNo.F1Target(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
% boxchart(1.6*ones(900,1),ThreshNo.precisionTarget(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
% boxchart(3.2*ones(900,1),ThreshNo.recallTarget(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
% boxchart(6*ones(900,1),ThreshNo.F1Interf(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
% boxchart(7.6*ones(900,1),ThreshNo.precisionInterf(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
% boxchart(9.2*ones(900,1),ThreshNo.recallInterf(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
% yyaxis right
% ylabel('\DeltaTIR (dB)')
% boxchart(13*ones(900,1),ThreshNo.deltaSNR_H(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
% 
% xticks([0, 1.6, 3.2, 6, 7.6, 9.2, 13])
% xticklabels({'F1 target','Precision target', 'Recall target', ...
%     'F1 interf','Precision interf', 'Recall interf','\DeltaTIR'})
%

%% TIR improvement

figure();
subplot(1,2,1)
datasetHistogram(ThreshYes.deltaSNR_H,-5:1:10,[-5 10],[0 0.4])
xlabel('\DeltaTIR (dB)')
%     title('\DeltaTIR');
set(gca,'FontSize',11);

subplot(1,2,2)
datasetHeatmap(ThreshYes.deltaSNR_H,'\DeltaTIR (dB)',[-5 10]);
%     title('\DeltaTIR');
set(gca,'FontSize',11);

figure();
subplot(1,2,1)
datasetHistogram(ThreshNo.deltaSNR_H,-5:1:10,[-5 10],[0 0.4])
xlabel('\DeltaTIR (dB)')
%     title('\DeltaTIR');
set(gca,'FontSize',11);

subplot(1,2,2)
datasetHeatmap(ThreshNo.deltaSNR_H,'\DeltaTIR (dB)',[-5 10]);
%     title('\DeltaTIR');
set(gca,'FontSize',11);

% Comparison boxplot
figure;
hold on
boxchart(-0.3*ones(900,1),ThreshYes.deltaSNR_H(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(0.3*ones(900,1),ThreshNo.deltaSNR_H(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
legend('with thresholding','without thresholding')
xticks([0])
xticklabels({''})
ylabel('\DeltaTIR (dB)')


%% TF plots

centerFreqsHz = [73.2236407563322	107.651956490872	...
    146.004400858780	188.728244811384	236.321739102815	...
    289.339924951948	348.401107001765	414.194064064606	...
    487.486081745890	569.131900623784	660.083684339571	...
    761.402123847844	874.268807325650	1000.00000000000	...
    1140.06199459580	1296.08821142304	1469.89824752649	...
    1663.51909705164	1879.20879030116	2119.48272716367	...
    2387.14301201846	2685.31113222345	3017.46436128696	...
    3387.47631126158	3799.66210728721	4258.82871111428	...
    4770.33098048645	5340.13411815414	5974.88323880703	...
    6681.98086522375];


% deltaTir_TF = TF_TIR_improvement(inputTargetSignal{1}, inputInterfSignal{1}, ...
%   outputTargetSignal_H{1}, outputInterfSignal_H{1}, AlgorithmStates, samplingRateHz, centerFreqsHz);
% figure;
% plotMonoSpectrogram(inputMixedSignal{1}, samplingRateHz, round(0.023 * samplingRateHz))

%% Control data set

ControlThreshNo = load('2022-06-05_02-48_meta_eval_data.mat');
displayMeanStd(ControlThreshNo.precisionTarget, ControlThreshNo.recallTarget, ...
    ControlThreshNo.precisionInterf, ControlThreshNo.recallInterf, ...
    ControlThreshNo.F1Target, ControlThreshNo.F1Interf, ControlThreshNo.deltaSNR_H)
ControlRandom = load('2022-06-05_11-33_meta_eval_data.mat');
displayMeanStd(ControlRandom.precisionTarget, ControlRandom.recallTarget,...
    ControlRandom.precisionInterf, ControlRandom.recallInterf, ...
    ControlRandom.F1Target, ControlRandom.F1Interf, ControlRandom.deltaSNR_H)
ControlEnhance = load('2022-06-06_03-11_meta_eval_data.mat');
displayMeanStd(ControlEnhance.precisionTarget, ControlEnhance.recallTarget,...
    ControlEnhance.precisionInterf, ControlEnhance.recallInterf, ...
    ControlEnhance.F1Target, ControlEnhance.F1Interf, ControlEnhance.deltaSNR_H)
ControlCancel = load('2022-06-06_11-04_meta_eval_data.mat');
displayMeanStd(ControlCancel.precisionTarget, ControlCancel.recallTarget,...
    ControlCancel.precisionInterf, ControlCancel.recallInterf, ...
    ControlCancel.F1Target, ControlCancel.F1Interf, ControlCancel.deltaSNR_H)

controlIndices = reshape(1:900,45,20);
controlIndices = controlIndices(1:2:45,:);
controlIndices = controlIndices(:,[2,3,9,11,13,15]);
controlIndices = controlIndices(:);

%% Symmetry

figure;
hold on
boxchart(1.7*ones(138,1),ThreshNo.precisionTarget(controlIndices),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(2.3*ones(138,1),ControlThreshNo.precisionTarget(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(3.7*ones(138,1),ThreshNo.recallTarget(controlIndices),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(4.3*ones(138,1),ControlThreshNo.recallTarget(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(7.7*ones(138,1),ThreshNo.precisionInterf(controlIndices),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(8.3*ones(138,1),ControlThreshNo.precisionInterf(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(9.7*ones(138,1),ThreshNo.recallInterf(controlIndices),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(10.3*ones(138,1),ControlThreshNo.recallInterf(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
legend('test set','control set')
xticks([2, 4, 8, 10])
xticklabels({'Target precision', 'Target recall', 'Interferer precision', 'Interferer recall'})
ylabel('Rate (1)')

figure;
hold on
boxchart(-0.3*ones(138,1),ThreshNo.deltaSNR_H(controlIndices),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(0.3*ones(138,1),ControlThreshNo.deltaSNR_H(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
legend('test set','control set')
xticks([0])
xticklabels({''})
ylabel('\DeltaTIR (dB)')

%% Random period

figure;
hold on
boxchart(1.7*ones(138,1),ControlThreshNo.precisionTarget(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(2.3*ones(138,1),ControlRandom.precisionTarget(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(3.7*ones(138,1),ControlThreshNo.recallTarget(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(4.3*ones(138,1),ControlRandom.recallTarget(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(7.7*ones(138,1),ControlThreshNo.precisionInterf(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(8.3*ones(138,1),ControlRandom.precisionInterf(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(9.7*ones(138,1),ControlThreshNo.recallInterf(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(10.3*ones(138,1),ControlRandom.recallInterf(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
legend('p_0 estimation','random p_0 assignment')
xticks([2, 4, 8, 10])
xticklabels({'Target precision', 'Target recall', 'Interferer precision', 'Interferer recall'})
ylabel('Rate (1)')

figure;
hold on
boxchart(-0.3*ones(138,1),ControlThreshNo.deltaSNR_H(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(0.3*ones(138,1),ControlRandom.deltaSNR_H(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
legend('p_0 estimation','random p_0 assignment')
xticks([0])
xticklabels({''})
ylabel('\DeltaTIR (dB)')

%% Enhancement, cancellation

figure;
hold on
boxchart(-0.8*ones(138,1),ControlThreshNo.deltaSNR_H(:),'BoxFaceColor','#0072BD','MarkerColor','#0072BD')
boxchart(0*ones(138,1),ControlEnhance.deltaSNR_H(:),'BoxFaceColor','#D95319','MarkerColor','#D95319')
boxchart(0.8*ones(138,1),ControlCancel.deltaSNR_H(:),'BoxFaceColor','#EDB120','MarkerColor','#EDB120')
legend('both','just enhancement', 'just cancellation')
xticks([0])
xticklabels({''})
ylabel('\DeltaTIR (dB)')

%%
function [m,s] = meanAndStd(input)
    fprintf('%.3f (%.3f)\n',mean(input,'all'),std(input,[],'all'));
end

function [] = displayMeanStd(precisionTarget, recallTarget, ...
    precisionInterf, recallInterf, F1Target, F1Interf, deltaSNR_H)

    fprintf('target precision - mean (std):\n')
    meanAndStd(precisionTarget)
    
    fprintf('target recall - mean (std):\n')
    meanAndStd(recallTarget)
    
    fprintf('interf precision - mean (std):\n')
    meanAndStd(precisionInterf)
    
    fprintf('interf recall - mean (std):\n')
    meanAndStd(recallInterf)
    
    fprintf('target F1 - mean (std):\n')
    meanAndStd(F1Target)
    
    fprintf('interf F1 - mean (std):\n')
    meanAndStd(F1Interf)
    
    fprintf('Delta TIR - mean (std):\n')
    meanAndStd(deltaSNR_H)

end