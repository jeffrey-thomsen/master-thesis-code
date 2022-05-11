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
function [processedSignal, AlgorithmStates, SimulationData] = ...
  blockFeedingRoutine(testSignal, BlockFeedingParameters, ...
  AlgorithmParameters, AlgorithmStates)

    warning('Blockwise storing of simulation data not compatible.')

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
    processedSignal = zeros(nBlocks*blockLength, size(testSignal,2));
    signalIndices = zeros(nBlocks,blockLength);
    for iBlock=1:nBlocks
        signalIndices(iBlock,:) = round((1:blockLength)+(iBlock-1)*blockLength);
        [processedSignal(signalIndices(iBlock,:),:), AlgorithmStates, ...
            SimulationData.Data(iBlock)] = ...
            speechEnhancement(testSignal(signalIndices(iBlock,:),:), ...
                AlgorithmParameters, AlgorithmStates);
    end
    SimulationData.blockLength = blockLength;
    SimulationData.nBlocks = nBlocks;
    SimulationData.signalIndices = signalIndices;
end