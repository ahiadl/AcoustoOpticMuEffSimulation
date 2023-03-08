close all
clear all %#ok<CLALL> 
% clc

%%
pd = PhantomDesigner();

uVars.wl   = 785; %[nm]
uVars.Vph  = 200;
uVars.musP = 0.4*ones(1,5); %[mm^-1]
uVars.mua  = 0.00258 * (2.^[1:5]); %[mm^-1]
uVars.numLayers = 1;
uVars.numPh = 5;
uVars.inkType = 'RoyalTalens';
uVars.layerThick = 60; %[mm]

pdRes = pd.createPhantomSet(uVars);

%%
aoSim = aoFluenceSim('MCX');

%% Simulate the fluence
uVarsAOSim = aoFluenceSim.createUserVars();

% General:
uVarsAOSim.savePath = './MuEffSimulation/Results';
uVarsAOSim.saveFlag = true;
uVarsAOSim.simName  = 'GradedAbsorption - measuredMuS';

%Graded Absorption
% uVarsAOSim.muaVec = [0.00516, 0.01032, 0.02064, 0.04128, 0.08256];
% uVarsAOSim.musVec = [0.89904, 0.89832, 0.89689, 0.89402, 0.88829];
% uVarsAOSim.musVec = [0.6737,  0.6737,  0.6737,  0.6737,  0.6737]; % measured with collimated transmission

uVarsAOSim.muaVec = pdRes.op.mua;
uVarsAOSim.musVec = pdRes.op.mus;

uVarsAOSim.trEnvSize = 1; %idxs

uVarsAOSim.g   = 0.555;%
uVarsAOSim.ref = 1.34;%

% FLuence:
uVarsAOSim.fluence.loadPhi       = false;
uVarsAOSim.fluence.phiPath       = [];
uVarsAOSim.fluence.srcPath       = "31-Jan-2023 08-59-25-Transmission-Uniform-NoSpacer.mat";
uVarsAOSim.fluence.detPath       = "31-Jan-2023 11-57-49-Reflection-Uniform-NoSpacer.mat";
uVarsAOSim.fluence.LightGeometry = 'Measured';
uVarsAOSim.fluence.srcParam      = [];
uVarsAOSim.fluence.nphoton       = 1e9; %1e8
uVarsAOSim.fluence.meshSize      = [60,60,60]; % lenX, lenY, lenZ
uVarsAOSim.fluence.meshRes       = 0.25; % [mm]

aoSim.simFluenceSet(uVarsAOSim);

% NOTICE: unfortunately, due to a bug in the mcx code, it is impossible to
% allocate gpuArray after using the mcx simulator. Therefore, we first
% simulate and then restarting the matlab and loading the mcx results.

%% Load the fluence and simulate AO:
uVarsAOSim = aoFluenceSim.createUserVars();

% General:
uVarsAOSim.savePath = './MuEffSimulation/Results/Graded Absorption';
uVarsAOSim.saveFlag = false;
uVarsAOSim.simName  = 'GradedAbsorptionFocused';

% FLuence:
uVarsAOSim.fluence.phiPath = './MuEffSimulation/Results/26-Feb-2023 22-28-37-GradedAbsorption-NewSet.mat';
uVarsAOSim.trEnvSize       = 15; %idxs

% VAOS:
uVarsAOSim.vaos.loadVaos       = false;
uVarsAOSim.vaos.vaosPath       = '';
uVarsAOSim.vaos.usPath         = "Analysis - FocusedAOITransducer-1.25MHz.mat";
% uVarsAOSim.vaos.usPath         = "Analysis - UnFocusedAOITransducer-1.25MHz";
uVarsAOSim.vaos.usFocalDist    = 60; %[mm] distance of the US focus from illumination plane
uVarsAOSim.vaos.N              = 251;
uVarsAOSim.vaos.fUS            = 1.25e6;
uVarsAOSim.vaos.pulseType      = 'measured';
uVarsAOSim.vaos.spacerLen      = 6.4; %mm
uVarsAOSim.vaos.spacerMaterial = 'PDMS';
uVarsAOSim.vaos.mathFluence    = false; 
uVarsAOSim.vaos.simSpeckle     = false;

uVarsAOSim.vaos.speckle.framesPerSig = 1000;
uVarsAOSim.vaos.speckle.sqncPerFrame = 10;
uVarsAOSim.vaos.speckle.batchSize    = 5;
uVarsAOSim.vaos.speckle.lambda       = 785e-9;
uVarsAOSim.vaos.speckle.gamma        = 1e-4;
uVarsAOSim.vaos.speckle.numOfGrain   = 1000;
uVarsAOSim.vaos.speckle.SBR          = 1e6;

% Meas Loader:
uVarsAOSim.ml.loadMeas     = true;
% uVarsAOSim.ml.measPath     = "D:/MuEff/OnlyAgar/Focused";
% uVarsAOSim.ml.measPath     = 'D:/MuEff/Uniform/GradedScatteringSet/Focused';
% uVarsAOSim.ml.measPath     = 'D:/MuEff/Uniform/GradedScatteringSet/UnFocused';
% uVarsAOSim.ml.measPath     = 'D:/MuEff/GradedAbsorption/Focused/WithPDMS';
% uVarsAOSim.ml.measPath     = 'D:/MuEff/GradedAbsorption/UnFocused/WithPDMS';
% uVarsAOSim.ml.measPath     = 'D:/MuEff/GradedAbsorption/Focused/WithAgar-SingleSet';
uVarsAOSim.ml.measPath     = './Measurements/MuEff/WithAgar-SingleSet';
uVarsAOSim.ml.measNameTemp = "Phantom-%d";
uVarsAOSim.ml.loadNew      = true;
uVarsAOSim.ml.idxLow       = 1; 
uVarsAOSim.ml.idxHigh      = 236;
uVarsAOSim.ml.phiHighToLow = true;
uVarsAOSim.ml.c            = [];
uVarsAOSim.ml.noiseIdxs    = 1:30;
uVarsAOSim.ml.alignToSig   = true;

% Graded Scattering
% uVarsAOSim.muEffIdx.log.sim   = 1345:1950;
uVarsAOSim.muEffIdx.log.speckle = 1:20;
% uVarsAOSim.muEffIdx.log.meas  = 1027:1340;

% Graded Absorption
uVarsAOSim.muEffLevels.sim     = [-4, -1];
uVarsAOSim.muEffLevels.speckle = 1:20;
uVarsAOSim.muEffLevels.meas    = [-1.99, -0.4];

aoSim.setVars(uVarsAOSim);
aoSim.config();
aoSim.extractMuEffFluence();
%%
% clc
close all
res = aoSim.simAndAnalyse();
%% Drafts:
clc;
close all
% aoSim.extractFluenceMuEff();
% aoSim.analyseMeas();
aoSim.compareMuEff();
aoSim.displayMuEffComparison();
aoSim.displayResults();
aoSim.displayResultsForIndexing();
aoSim.ml.displayResults
