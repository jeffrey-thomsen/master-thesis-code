function testedInputSignal = testInputSignal(inputSignal)    
    % run some basic checks on the input signal
    assert(ismatrix(inputSignal), "input signal must be 2D matrix")
    assert(isreal(inputSignal), "input signal must be real-valued")
    assert(all(abs(inputSignal(:))<=1), ...
        "input signal must be normalized to 1/-1")
    
    % make sure it has the right dimensions (n-by-2)
    if (size(inputSignal,2)~=2)
       inputSignal = transpose(inputSignal);
    end
    
    testedInputSignal = inputSignal;
end