function testSignal = testSignalGenerator(TestSignalParameters)
    
    if nargin>0
        nSamples = TestSignalParameters.nSamples;
    else
        nSamples = 1000;
    end

    testSignal = rand(nSamples,2)-0.5; % default placeholder for now

end