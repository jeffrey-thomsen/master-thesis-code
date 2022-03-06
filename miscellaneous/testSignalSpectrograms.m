% small script to visualize the speech signals I want to use for evaluating
% my speech enhancement algorithm

clear
close all

fileNames{1}='sp01.wav';
fileNames{2}='sp30.wav';
fileNames{3}='intelligence_16.wav';
fileNames{4}='storylines_16.wav';
fileNames{5}='p298_097.wav';
fileNames{6}='p313_256.wav';
fileNames{7}='p360_253.wav';
fileNames{8}='2078-142845-0002.flac';


for iFile = 1:length(fileNames)
    [signal{iFile},fs{iFile}] = ...
        audioread(fileNames{iFile});
end

for iFile = 1:length(fileNames)
    figure;
    spectrogram(signal{iFile},hamming(0.2*fs{iFile}),[],[],fs{iFile});
    xlim([0.1, 0.5]);
    title(['signal: ', fileNames{iFile}, ' - spectrogram'])
end