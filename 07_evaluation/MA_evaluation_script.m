
tic

dateString = '2022-05-27_02-55';

load([dateString,'_evaluation_data.mat']);
load([dateString,'_signal_data.mat'])
load([dateString,'_simulation_data.mat'])

centerFreqsHz = AlgorithmStates.L.GammatoneStates.analyzer.center_frequencies_hz;

doaHistogram_subplots;

[precisionTarget, precisionInterf, recallTarget, recallInterf, ...
    F1Target, F1Interf] = evaluationPlots(speakerCombinations, ...
    anglePermutations, precision, recall, deltaSNR_H);

deltaSnr_TF = TF_SNR_improvement(inputTargetSignal{1}, inputInterfSignal{1}, ...
  outputTargetSignal_H{1}, outputInterfSignal_H{1}, AlgorithmStates, samplingRateHz, centerFreqsHz);
figure;
plotMonoSpectrogram(inputMixedSignal{1}, samplingRateHz, round(0.023 * samplingRateHz))


fprintf('target precision:')
meanAndStd(precisionTarget)

fprintf('target recall:')
meanAndStd(recallTarget)

fprintf('interf precision:')
meanAndStd(precisionInterf)

fprintf('interf recall:')
meanAndStd(recallInterf)

fprintf('Delta SNR:')
meanAndStd(deltaSNR_H)

toc

function [m,s] = meanAndStd(input)
    mean(input,'all')
    std(input,[],'all')
end