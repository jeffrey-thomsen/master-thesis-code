        % compute IBM and extract glimpse masks from simulation data
        [inputTargetSubbandSignals, ~] = subbandDecompositionBinaural(...
            inputTargetSignal{iSp,jAn}, AlgorithmStates);
        [inputInterfSubbandSignals, ~] = subbandDecompositionBinaural(...
            inputInterfSignal{iSp,jAn}, AlgorithmStates);
        
        inputSnr{iSp,jAn}.L = computeTimeFreqSnr(inputTargetSubbandSignals.L, ...
            inputInterfSubbandSignals.L, samplingRateHz, nBands);
        inputSnr{iSp,jAn}.R = computeTimeFreqSnr(inputTargetSubbandSignals.R, ...
        inputInterfSubbandSignals.R, samplingRateHz, nBands);        
        
        % compute TF gray mask SNR improvement
        [outputTargetSubbandSignals_H, ~] = subbandDecompositionBinaural(...
            outputTargetSignal_H{iSp,jAn}, AlgorithmStates);
        [outputInterfSubbandSignals_H, ~] = subbandDecompositionBinaural(...
            outputInterfSignal_H{iSp,jAn}, AlgorithmStates);

        outputSnr{iSp,jAn}.L = computeTimeFreqSnr(outputTargetSubbandSignals_H.L, ...
            outputInterfSubbandSignals_H.L, samplingRateHz, nBands);
        outputSnr{iSp,jAn}.R = computeTimeFreqSnr(outputTargetSubbandSignals_H.R, ...
            outputInterfSubbandSignals_H.R, samplingRateHz, nBands);