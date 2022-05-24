% A function to emulate the Bruemann2018 algorithm, as a means of
% comparision within my framework. This will only work with target angle 0°
function azimuthDegCells = bruemannAzimuthEstimation(subbandSignalArray, ...
  AlgorithmParameters, AlgorithmStates)

    warning('Bruemann simple azimiuth estimation is not implemented for block-wise processing, only whole signal!')
    
    nBands = AlgorithmParameters.Gammatone.nBands;
    azimuthDegCells = cell(1, nBands);

    centerFreqsHz = AlgorithmStates.L.GammatoneStates.analyzer.center_frequencies_hz;

%     angleLim = 5;
%     timeLim = 0.7e-3*sin((angleLim*2*pi)/360);

    for iBand = 1:nBands
        subbandSignal.L = subbandSignalArray.L(:,iBand);
        subbandSignal.R = subbandSignalArray.R(:,iBand);
        
        ITF = subbandSignal.L./subbandSignal.R;
        % ITF = R.*conj(L);
        % atan(imag(ITF)./real(ITF))
        LR = bsxfun(@rdivide, angle(ITF), (2 * pi * centerFreqsHz(iBand)).');
        [ITD, ~] = firstOrderLowPass(LR, 0, AlgorithmParameters);
        assert(max(abs(ITD))<=1, 'ITD values lie outside of maximum range')
        azimuthDeg = (360 / (2 * pi)) * asin(ITD / 0.7e-3);
        azimuthDegCells{iBand} = azimuthDeg;
    end

end