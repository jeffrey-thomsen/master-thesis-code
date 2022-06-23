# Speech Enhancement for Hearing Aids Using a Combination of Binaural and Periodicity Features

This is the MATLAB code of the speech enhancement algorithm developed during my master thesis at Universität Oldenburg in 2022.

To try it out, simply clone or download the repository, and run the *thomsen2022_main.m* script.

Notes:
- You must add the folder of the repository and all its subfolders to your MATLAB search path.
- You need to install the [SOFA API for Matlab/Octave](https://sourceforge.net/projects/sofacoustics/).

The main function of the speech enhancement algorithm is *01_base/speechEnhancement.m*. Its subfunctions are stored in the folders *03_filterbank*, *04_periodicityAnalysis*, *05_directionOfArrivalEstimation* and *06_enhancementStage*.

The folder *02_testSignalGeneration* contains a few test stimuli and functions used by *thomsen2022_main.m* for generating binaural multi-talker test signals.

The folder *Hohmann2002* contains the original implementation of the gammatone filterbank that is used by the *subbandDecomposition.m* and *subbandResynthesis.m* functions. It is taken from the [Auditory Modeling Toolbox](https://sourceforge.net/projects/amtoolbox/).

*testSuite.m* and *TestAlgorithm.m* contain unit tests that were used during the development of the algorithm.
