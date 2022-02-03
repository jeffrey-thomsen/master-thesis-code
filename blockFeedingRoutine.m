function processedSignal = blockFeedingRoutine(testSignal, ...
  BlockFeedingParameters, AlgorithmParameters)

    if nargin > 1
        blockLength = BlockFeedingParameters.blockLength;
    else
        blockLength = 100;
    end
    assert(blockLength<length(testSignal))
    nBlocks = floor(length(testSignal)/blockLength);

    
    % Cut signal into specified block sizes and run speech enhancement
    % algorithm block-by-block
    processedSignal = zeros(size(testSignal));
    for iBlock=1:nBlocks
        signalIndices = round((1:blockLength)+(iBlock-1)*blockLength);
        processedSignal(signalIndices,:) = ...
            speechEnhancement(testSignal(signalIndices,:));
    end
end