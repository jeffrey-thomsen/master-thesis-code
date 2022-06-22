% create a subset of data from the main set that can be compared to the control set
% [90 0 -30]
a.AlgorithmParameters = AlgorithmParameters;
a.anglePermutations = anglePermutations(controlAngles,:);
a.deltaSNR_H = deltaSNR_H(controlSpeakers,controlAngles);
a.deltaSNR_S = deltaSNR_S(controlSpeakers,controlAngles);
a.maskAppliedTarget = maskAppliedTarget(controlSpeakers,controlAngles);
a.maskAppliedInterf = maskAppliedInterf(controlSpeakers,controlAngles);
a.MetaData = MetaData;
a.precision = precision(controlSpeakers,controlAngles);
a.recall = recall(controlSpeakers,controlAngles);
a.speakerCombinations = speakerCombinations(controlSpeakers,:);
a.TestSignalParameters = TestSignalParameters;
