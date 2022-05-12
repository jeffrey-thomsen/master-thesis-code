% Compute SNR improvement yielded from speech enhancement algorithm by
% computing and then subtracting the SNR of input and output target and
% interferer signals
function snrImprovement = computeSnrImprovement(targetSignal, interfSignal, ...
    enhancedSignalTarget, enhancedSignalInterf)

    snrInit = 10 * log10(sum(abs(targetSignal(:)).^2) ./ ...
                         sum(abs(interfSignal(:)).^2));
    snrPost = 10 * log10(sum(abs(enhancedSignalTarget(:)).^2) ./ ...
                         sum(abs(enhancedSignalInterf(:)).^2));

    snrImprovement = snrPost - snrInit;

end