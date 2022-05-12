% Compute precision and recall values by comparing logical matrices,
% denoting time-frequency glimpses chosen as target/interferer by the
% algorithm to logical matrices of the ideal binary mask (IBM), denoting 
% the theoretical optimum performance of the estimation stage of the
% algorithm
function [precision, recall] = compareToIbm(...
  maskTarget, maskInterf, ibmTarget, ibmInterf)

    nTruePositivesTarget = nnz(maskTarget.L & ibmTarget.L) +...
                           nnz(maskTarget.R & ibmTarget.R);
    nTruePositivesInterf = nnz(maskInterf.L & ibmInterf.L) +...
                           nnz(maskInterf.R & ibmInterf.R);

    nFalsePositivesTarget = nnz(maskTarget.L & ~ibmTarget.L) +...
                            nnz(maskTarget.R & ~ibmTarget.R);
    nFalsePositivesInterf = nnz(maskInterf.L & ~ibmInterf.L) +...
                            nnz(maskInterf.R & ~ibmInterf.R);

    nFalseNegativesTarget = nnz(~maskTarget.L & ibmTarget.L) +...
                            nnz(~maskTarget.R & ibmTarget.R);
    nFalseNegativesInterf = nnz(~maskTarget.L & ibmTarget.L) +...
                            nnz(~maskTarget.R & ibmTarget.R);

    precision.Target = nTruePositivesTarget / ...
        (nTruePositivesTarget + nFalsePositivesTarget);
    precision.Interf = nTruePositivesInterf / ...
        (nTruePositivesInterf + nFalsePositivesInterf);
    recall.Target = nTruePositivesTarget / ...
        (nTruePositivesTarget + nFalseNegativesTarget);
    recall.Interf = nTruePositivesInterf / ...
        (nTruePositivesInterf + nFalseNegativesInterf);
    
    precision.Total = (nTruePositivesTarget + nTruePositivesInterf) / ...
        (nTruePositivesTarget + nFalsePositivesTarget + ...
         nTruePositivesInterf + nFalsePositivesInterf);
end