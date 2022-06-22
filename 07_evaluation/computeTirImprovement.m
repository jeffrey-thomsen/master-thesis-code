% Compute TIR improvement yielded from speech enhancement algorithm by
% computing and then subtracting the TIR of input and output target and
% interferer signals
function tirImprovement = computeTirImprovement(targetSignal, interfSignal, ...
    enhancedSignalTarget, enhancedSignalInterf)

    tirIn  = 10 * log10(sum(abs(targetSignal(:)).^2) ./ ...
                         sum(abs(interfSignal(:)).^2));
    tirOut = 10 * log10(sum(abs(enhancedSignalTarget(:)).^2) ./ ...
                         sum(abs(enhancedSignalInterf(:)).^2));

    tirImprovement = tirOut - tirIn;

end