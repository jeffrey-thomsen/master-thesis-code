function [testSignal, targetSignal, interfSignal, fsHrtf, ...
  anglePermutations, speakerCombinations] = ...
  testSignalGenerator(TestSignalParameters, hrtf)

% check desired test signal type (keep compatible with old unit tests)
if nargin == 0
    TestSignalParameters.testSignalType = 'WhiteNoise';
elseif nargin == 1
   if ~isfield(TestSignalParameters, 'testSignalType')
       TestSignalParameters.testSignalType = 'WhiteNoise';
   end
   if ~(strcmp(TestSignalParameters.testSignalType, 'WhiteNoise'))
       hrtf = SOFAload('HRIR_KEMAR_DV0001_4.sofa',[5 2],'R');
   end
end

if ~(strcmp(TestSignalParameters.testSignalType, 'WhiteNoise'))
    fsHrtf = hrtf.Data.SamplingRate;
end

switch TestSignalParameters.testSignalType
    case 'WhiteNoise'

        if isfield(TestSignalParameters, 'nSamples')
            nSamples = TestSignalParameters.nSamples;
        else
            nSamples = 32000;
        end
    
        testSignal = rand(nSamples,2)-0.5; % default placeholder for now
    
    case 'DAGA'       
    
        % load clean speech signals
        [targetSignal, fsTarget] = audioread('pathological_16.wav');
        [interfSignal, fsInterf] = audioread('peaches_16.wav');
    
        % adjust sampling rate to HRTF
        targetSignal = resample(targetSignal, fsHrtf, fsTarget);
        interfSignal = resample(interfSignal, fsHrtf, fsInterf);
    
        % equalize lengths
        if length(targetSignal) > length(interfSignal)
            targetSignal = targetSignal(1:length(interfSignal));
        elseif length(targetSignal) < length(interfSignal)
            interfSignal = interfSignal(1:length(targetSignal));
        end
    
        % convolve with HRTF at different incidence angles
        targetSignal = SOFAspat(targetSignal, hrtf, 0, 0);
        interfSignal = SOFAspat(interfSignal, hrtf, 60, 0);
        
        % equalize levels - Target-to-interferer energy ratio 0dB
        targetSignal = targetSignal / std(targetSignal(:));
        interfSignal = interfSignal / std(interfSignal(:));
        
        % add and normalize test signal
        mixedSignal = targetSignal + interfSignal;
        
        scalingFactor = max(max(mixedSignal));
        
        testSignal  = mixedSignal ./scalingFactor;
        targetSignal = targetSignal./scalingFactor;
        interfSignal = interfSignal./scalingFactor;

    case 'Battery'

        % read parameters
        targetAngles = TestSignalParameters.targetAngles;
        nAngles = numel(targetAngles);
        speakerIdVector = TestSignalParameters.speakerIds;
        nSpeakers = TestSignalParameters.nSpeakers;
        assert(nAngles>=nSpeakers, 'Number of angles must be at least the number of speakers');
%         snr = TestSignalParameters.snr;
    
        % generate battery of possible voice combinations
    
        % load all speech signals
        iSource = 0;
        for speakerId = speakerIdVector
            iSource = iSource+1;
            [speech.sig, speech.fs] = ...
                audioread(strcat(speakerId, '.flac'));
            speechCell{iSource} = speech;
            speechLength(iSource) = length(speech.sig);
        end
        nSources = iSource;

        % equalize lengths
        [~, minInd] = min(speechLength);
        minLength = speechLength(minInd);
        for iSource = 1:nSources
            speechCell{iSource}.sig = speechCell{iSource}.sig(1:minLength);
        end
        
        % resample signals to HRTF sampling rate
        speechCellResampled = cell([nSources, 1]);
        for iSource = 1:nSources
            speechCellResampled{iSource} = ...
                resample(speechCell{iSource}.sig, fsHrtf, ...
                speechCell{iSource}.fs);
        end
        
        % spatialize each speech signal with every desired DOA
        speechCellSpatialized = cell([nSources, nAngles]);
        for iSource = 1:nSources
            for jAngle = 1:nAngles
                speechCellSpatialized{iSource, jAngle} = ...
                    SOFAspat(speechCellResampled{iSource}, hrtf, ...
                    targetAngles(jAngle), 0);
                % equalize levels - Target-to-interferer energy ratio 0dB
                speechCellSpatialized{iSource, jAngle} = ...
                    speechCellSpatialized{iSource, jAngle} / ...
                    std(speechCellSpatialized{iSource, jAngle}(:));
            end
        end
        
        % MISSING: VARIABLE SNR

        % compute all possible combinations of speakers and angles
        nk = nchoosek(1:nAngles, nSpeakers);
        p = zeros(0, nSpeakers);
        for i = 1:size(nk, 1)
            pi = perms(nk(i, :));
            p = unique([p; pi], 'rows');
        end
        anglePermutations = p;%perms(1:nAngles);
        nAnglePerms = size(anglePermutations, 1);
        speakerCombinations = nchoosek(1:nSources, nSpeakers);
        nSpeakerCombos = size(speakerCombinations, 1);
        
        % add up and normalize test signal
        mixedSignalCell  = cell([nSpeakerCombos, nAnglePerms]);
        targetSignalCell = cell([nSpeakerCombos, nAnglePerms]);
        interfSignalCell = cell([nSpeakerCombos, nAnglePerms]);
        for iSpeakerCombo = 1:nSpeakerCombos
            for jAnglePerm = 1:nAnglePerms
                targetSignalCell{iSpeakerCombo, jAnglePerm} = ...
                        speechCellSpatialized{...
                        speakerCombinations(iSpeakerCombo, 1), ...
                        anglePermutations(jAnglePerm, 1)...
                        };
                interfSignalCell{iSpeakerCombo, jAnglePerm} = ...
                        speechCellSpatialized{...
                        speakerCombinations(iSpeakerCombo, 2), ...
                        anglePermutations(jAnglePerm, 2)...
                        };
                for kSpeaker = 3:nSpeakers % if >2 speakers present
                    interfSignalCell{iSpeakerCombo, jAnglePerm} = ...
                        interfSignalCell{iSpeakerCombo, jAnglePerm} + ...
                        speechCellSpatialized{...
                        speakerCombinations(iSpeakerCombo, kSpeaker), ...
                        anglePermutations(jAnglePerm, kSpeaker)...
                        };
                end

                mixedSignalCell{iSpeakerCombo, jAnglePerm} = ...
                    targetSignalCell{iSpeakerCombo, jAnglePerm} + ...
                    interfSignalCell{iSpeakerCombo, jAnglePerm};

                scalingFactor = max(max(mixedSignalCell{iSpeakerCombo, jAnglePerm}));
                mixedSignalCell{iSpeakerCombo, jAnglePerm} = ...
                    mixedSignalCell{iSpeakerCombo, jAnglePerm} / scalingFactor;
                targetSignalCell{iSpeakerCombo, jAnglePerm} = ...
                    targetSignalCell{iSpeakerCombo, jAnglePerm} / scalingFactor;
                interfSignalCell{iSpeakerCombo, jAnglePerm} = ...
                    interfSignalCell{iSpeakerCombo, jAnglePerm} / scalingFactor;
            end
        end

        testSignal = mixedSignalCell;
        targetSignal = targetSignalCell;
        interfSignal = interfSignalCell;
        anglePermutations = targetAngles(anglePermutations);
        speakerCombinations = speakerIdVector(speakerCombinations);%nchoosek(speakerIdVector, nAngles);

        % MIGHT WANT: global scaling factor for all test signals

end


end