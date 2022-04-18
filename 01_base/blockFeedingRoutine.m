% Divide a binaural signal into blocks, to reduce memory load, and run the
% Thomsen2022 speech enhancement algorithm on each block.
% Block length can also be chosen to be 1 to showcase cauality of the
% algorithm.
% testSignal - Nx2 matrix containing real-valued binaural input signal
% processedSignal - Nx2 matrix with the real-valued processed signal
% ..Parameters - structs containing parametric information for the
% simulation
% AlgorithmStates - struct containing e.g. filter states that are updated
% every sample and need to be updated for each processed block
function [processedSignal, AlgorithmStates, SimulationData] = blockFeedingRoutine(testSignal, ...
  BlockFeedingParameters, AlgorithmParameters, AlgorithmStates)

    warning('Do not have propoer simulation data gathering in block feeding routing yet!')

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
        [processedSignal(signalIndices,:), AlgorithmStates, SimulationData] = ...
            speechEnhancement(testSignal(signalIndices,:), ...
                AlgorithmParameters, AlgorithmStates);
    end
end