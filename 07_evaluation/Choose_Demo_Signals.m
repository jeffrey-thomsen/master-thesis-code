%% Main set with thresh yes/no

% 2022-05-27_02-55 yes
% 2022-05-29_23-03 no

samplingRateHz = 13400;

iSp = 5;%25,40,34,20
jAn = 6;

scaleFactor = max(...
[max(abs(inputMixedSignal{iSp,jAn}),           [], 'all'),...
 max(abs(outputMixedSignalThreshYes{iSp,jAn}), [], 'all'),...
 max(abs(outputMixedSignalThreshNo{iSp,jAn}),  [], 'all')]);

box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
sound(inputMixedSignal{iSp,jAn}./scaleFactor, samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
sound(outputMixedSignalThreshYes{iSp,jAn}./scaleFactor, samplingRateHz);
box3 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box3);
sound(outputMixedSignalThreshNo{iSp,jAn}./scaleFactor, samplingRateHz);


%5/6, 5/10

audiowrite(['main_input','.wav'], inputMixedSignal{iSp,jAn}./scaleFactor, samplingRateHz);
audiowrite(['main_output_threshyes','.wav'], outputMixedSignalThreshYes{iSp,jAn}./scaleFactor, samplingRateHz);
audiowrite(['main_output_threshno','.wav'], outputMixedSignalThreshNo{iSp,jAn}./scaleFactor, samplingRateHz);

%% Control set with cancel, enhance, both

% 2022-06-05_02-48 both
% 2022-06-06_03-11 enhance
% 2022-06-06_11-04 cancel

iSp = 3;%12
jAn = 6;%4

scaleFactor = max(...
[max(abs(inputMixedSignalControl{iSp,jAn}),  [], 'all'),...
 max(abs(outputMixedSignalCancel{iSp,jAn}),  [], 'all'),...
 max(abs(outputMixedSignalEnhance{iSp,jAn}), [], 'all'),...
 max(abs(outputMixedSignalBoth{iSp,jAn}),    [], 'all')]);

box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
sound(inputMixedSignalControl{iSp,jAn}./scaleFactor, samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
sound(outputMixedSignalCancel{iSp,jAn}./scaleFactor, samplingRateHz);
box3 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box3);
sound(outputMixedSignalEnhance{iSp,jAn}./scaleFactor, samplingRateHz);
box4 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box4);
sound(outputMixedSignalBoth{iSp,jAn}./scaleFactor, samplingRateHz);

audiowrite(['control_input','.wav'], inputMixedSignalControl{iSp,jAn}./scaleFactor, samplingRateHz);
audiowrite(['control_output_cancel','.wav'], outputMixedSignalCancel{iSp,jAn}./scaleFactor, samplingRateHz);
audiowrite(['control_output_enhance','.wav'], outputMixedSignalEnhance{iSp,jAn}./scaleFactor, samplingRateHz);
audiowrite(['control_output_both','.wav'], outputMixedSignalBoth{iSp,jAn}./scaleFactor, samplingRateHz);