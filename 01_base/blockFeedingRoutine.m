% Divide a binaural signal into blocks reducing memory load and run the
% Thomsen2022 speech enhancement algorithm on each block.
% Block length can also be chosen to be 1 to showcase cauality of the
% algorithm.
% testSignal - Nx2 matrix containing real-valued binaural input signal
% processedSignal - Nx2 matrix with the real-valued processed signal
function [processedSignal, AlgorithmStates] = blockFeedingRoutine(testSignal, ...
  BlockFeedingParameters, AlgorithmParameters, AlgorithmStates)

    if nargin > 1 && isfield(BlockFeedingParameters,'blockLength')
        blockLength = BlockFeedingParameters.blockLength;
    else
        blockLength = 100;
    end
    assert(blockLength<length(testSignal), ...
        "Block length must be shorter than signal")
    nBlocks = floor(size(testSignal,1)/blockLength);

    % Cut signal into specified block sizes and run speech enhancement
    % algorithm block-by-block
    processedSignal = zeros(size(testSignal));
    for iBlock=1:nBlocks
        signalIndices = round((1:blockLength)+(iBlock-1)*blockLength);
        [processedSignal(signalIndices,:), AlgorithmStates] = ...
            speechEnhancement(testSignal(signalIndices,:), ...
                AlgorithmParameters, AlgorithmStates);
    end
end