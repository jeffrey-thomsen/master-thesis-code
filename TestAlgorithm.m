%% testing algorithm parameters
classdef TestAlgorithm < matlab.unittest.TestCase
    properties (TestParameter)
        ChenP0Detection = {true, false};
        coherenceMask   = {true, false};
        azimuthPooling  = {true, false};
        cfr0mask        = {true, false};
        DOAProcessing   = {true, false};
        Cancellation    = {true, false};
        Enhancement     = {true, false};
    end

    methods (Test)
        function testCancellation(testCase, ChenP0Detection, ...
          coherenceMask, azimuthPooling, cfr0mask, ...
          DOAProcessing, Cancellation, Enhancement)
        % test if the speech enhancement algorithm will run with all
        % combinations of processing parameters

            % Setup
            AlgorithmParameters = AlgorithmParametersConstructor();

            AlgorithmParameters.ChenP0Detection = ChenP0Detection;
            AlgorithmParameters.coherenceMask   = coherenceMask;
            AlgorithmParameters.azimuthPooling  = azimuthPooling;
            AlgorithmParameters.cfr0mask        = cfr0mask;
            AlgorithmParameters.DOAProcessing   = DOAProcessing;
            AlgorithmParameters.Cancellation    = Cancellation;
            AlgorithmParameters.Enhancement     = Enhancement;

            lookuptable = load('2022-03-07_itd_lookuptable_annotated.mat');
            lookuptable = lookuptable.lookuptable.lookuptable;
            AlgorithmParameters.lookuptable = lookuptable;
            [AlgorithmStates, AlgorithmParameters.Gammatone.nBands] = ...
                AlgorithmStatesConstructor(AlgorithmParameters);

            % create test signal
            input = rand(500,2)-0.5;
            input = testInputSignal(input);

            % Exercise
            [output, ~, simulationData] = ...
                speechEnhancement(input, AlgorithmParameters, AlgorithmStates);

            % Verify

            verifySize(testCase,output,size(input))
            
            if azimuthPooling && ~coherenceMask && (Cancellation || Enhancement)
                for iBand = 16:30
                    verifyEqual(testCase, ...
                        simulationData.azimuthDegCells{iBand}, ...
                        simulationData.azimuthDegCells{iBand-1})
                end
            end
        end
    end
end