% Main set with thresh yes/no

% 2022-05-27_02-55 yes
% 2022-05-29_23-03 no

iSp = 5;%25,40,34,20
jAn = 6;
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
soundsc(inputMixedSignal{iSp,jAn}, samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
soundsc(outputMixedSignalThreshYes{iSp,jAn}, samplingRateHz);
box3 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box3);
soundsc(outputMixedSignalThreshNo{iSp,jAn}, samplingRateHz);
%5/6, 5/10

% Control set with cancel, enhance, both

% 2022-06-05_02-48 both
% 2022-06-06_03-11 enhance
% 2022-06-06_11-04 cancel

iSp = 12;%12
jAn = 4;%4
box1 = msgbox('Play original signal (Ensure volume is adequately set)');
waitfor(box1);
soundsc(inputMixedSignalControl{iSp,jAn}, samplingRateHz);
box2 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box2);
soundsc(outputMixedSignalCancel{iSp,jAn}, samplingRateHz);
box3 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box3);
soundsc(outputMixedSignalEnhance{iSp,jAn}, samplingRateHz);
box4 = msgbox(['Play resynthesized signal. To replay, just rerun last' ...
    ' section of script (adjusting filter variables if necessary)']);
waitfor(box4);
soundsc(outputMixedSignalBoth{iSp,jAn}, samplingRateHz);