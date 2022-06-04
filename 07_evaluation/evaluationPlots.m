function [precisionTarget, precisionInterf, recallTarget, recallInterf, ...
    F1Target, F1Interf] = evaluationPlots(varargin)

    if nargin == 1
        EvalStruct = varargin{1};
    elseif nargin == 5
        EvalStruct.speakerCombinations = varargin{1};
        EvalStruct.anglePermutations = varargin{2};
        EvalStruct.precision = varargin{3};
        EvalStruct.recall = varargin{4};
        EvalStruct.deltaSNR_H = varargin{5};
    else
        error('Must specify 1 struct or 5 variables as input!')
    end
    %% Evaluation evaluation
    
    nSpeakerCombos = size(EvalStruct.speakerCombinations, 1);
    nAnglePerms = size(EvalStruct.anglePermutations, 1);
    
    
    figure;
    subplot(1,2,1)
    histogram(EvalStruct.deltaSNR_H);
    title('delta SNR');
    subplot(1,2,2)
    heatmap(EvalStruct.deltaSNR_H);
    title('delta SNR');
    
    
    for iSp=1:nSpeakerCombos
        for jAn=1:nAnglePerms
            precisionTarget(iSp,jAn)=EvalStruct.precision{iSp,jAn}.Target;
        end
    end
    figure;
    subplot(1,2,1)
    histogram(precisionTarget);
    title('target precision');
    subplot(1,2,2)
    heatmap(precisionTarget);
    title('target precision');
    for iSp=1:nSpeakerCombos
        for jAn=1:nAnglePerms
            precisionInterf(iSp,jAn)=EvalStruct.precision{iSp,jAn}.Interf;
        end
    end
    figure;
    subplot(1,2,1)
    histogram(precisionInterf);
    title('interf precision');
    subplot(1,2,2)
    heatmap(precisionInterf);
    title('interf precision');
    
    
    
    for iSp=1:nSpeakerCombos
        for jAn=1:nAnglePerms
            recallTarget(iSp,jAn)=EvalStruct.recall{iSp,jAn}.Target;
        end
    end
    figure;
    subplot(1,2,1)
    histogram(recallTarget);
    title('target recall');
    subplot(1,2,2)
    heatmap(recallTarget);
    title('target recall');
    for iSp=1:nSpeakerCombos
        for jAn=1:nAnglePerms
            recallInterf(iSp,jAn)=EvalStruct.recall{iSp,jAn}.Interf;
        end
    end
    figure;
    subplot(1,2,1)
    histogram(recallInterf);
    title('interf recall');
    subplot(1,2,2)
    heatmap(recallInterf);
    title('interf recall');
    
    %% F1 score
    
    F1Target = 2 .* (precisionTarget.*recallTarget)./(precisionTarget+recallTarget);
    figure;
    subplot(1,2,1)
    histogram(F1Target);
    title('F1 target');
    subplot(1,2,2)
    heatmap(F1Target);
    title('F1 target');
    
    F1Interf = 2 .* (precisionInterf.*recallInterf)./(precisionInterf+recallInterf);
    figure;
    subplot(1,2,1)
    histogram(F1Interf);
    title('F1 interf');
    subplot(1,2,2)
    heatmap(F1Interf);
    title('F1 interf');

end