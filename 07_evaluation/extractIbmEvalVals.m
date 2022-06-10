function [precisionTarget, precisionInterf, recallTarget, recallInterf, ...
    F1Target, F1Interf] = extractIbmEvalVals(varargin)

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
    
    for iSp=1:nSpeakerCombos
        for jAn=1:nAnglePerms
            precisionTarget(iSp,jAn)=EvalStruct.precision{iSp,jAn}.Target;
            precisionInterf(iSp,jAn)=EvalStruct.precision{iSp,jAn}.Interf;
            recallTarget(iSp,jAn)=EvalStruct.recall{iSp,jAn}.Target;
            recallInterf(iSp,jAn)=EvalStruct.recall{iSp,jAn}.Interf;
        end
    end

   
    
    %% F1 score
    
    F1Target = 2 .* (precisionTarget.*recallTarget)./(precisionTarget+recallTarget);
%     figure;
%     subplot(1,2,1)
%     histogram(F1Target);
%     title('F1 target');
%     subplot(1,2,2)
%     heatmap(F1Target);
%     title('F1 target');
    
    F1Interf = 2 .* (precisionInterf.*recallInterf)./(precisionInterf+recallInterf);
%     figure;
%     subplot(1,2,1)
%     histogram(F1Interf);
%     title('F1 interf');
%     subplot(1,2,2)
%     heatmap(F1Interf);
%     title('F1 interf');

end