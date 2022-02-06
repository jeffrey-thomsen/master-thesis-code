% Ensures that the signal used for testing the Thomsen2022 speech
% enhancement algorithm has the correct dimensions (Nx2) and transposes the
% matrix if necessary. Asserts that the signal has two channels, is
% real-valued and is normalized to 1/-1
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
    assert(size(inputSignal,2)==2, ...
        "input signal must contain exactly two channels")
    
    testedInputSignal = inputSignal;
end