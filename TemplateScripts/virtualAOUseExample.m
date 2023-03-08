close all;
clear all;
clc;

vAOSim = virtualAOSim();
usPath = ".\AcoustoOpticSystem\Measurements\Transducer Pressure Field\Calibrated\Focused\Analysis - FocusedAOITransducer-1.25MHz.mat";
vAOSim.loadUS(usPath);
%%
clc
close all

uVars = vAOSim.createUserVars();

uVars.debug = false;
uVars.displayDebug = false;
uVars.debugTime = false;

uVars.N   = 251;
uVars.fUS = 1.25e6;
uVars.pulseType = 'measured'; %delta, measured, 

uVars.spacerLen = 6.4; %mm
uVars.spacerMaterial = 'PDMS';

uVars.usDistFromInt = 30; %[mm] 

uVars.useCustomUSParams = false;

uVars.muEffVec      = [0.074; 0.106; 0.152; 0.219; 0.323];

vAOSim.setVars(uVars)
simVars = vAOSim.getVars();

vAOSim.calcSimDimAndSpace();
vAOSim.createMathematicalFluence();
vAOSim.alignAndInterpUS();
vAOSim.createAcousticInterface();
vAOSim.createPulses();
vAOSim.buildSPMatrix();
vAOSim.buildHadMatrix();
vAOSim.buildHadInvMat();
vAOSim.reconAll(vAOSim.phiMath);
vAOSim.displayResults();

resSpeckle = vAOSim.speckleSim(vAOSim.phiMath, true);
%% 
% clc
% close all
vAOSim.createMathematicalFluence();

vAOSim.reconAll(vAOSim.phiMath);
% resSpeckle = vAOSim.speckleSim(vAOSim.phiMath, true);

%% Load MCX Simulations:
resSim    = load("13-Nov-2022 14-22-49-FocusedTrans.mat");

%% Extract MCX Simulations:
envSize = 0;
trSize = 2*envSize+1;
phiSimRaw = zeros(5, 241, trSize, trSize);
varsSim = resSim.simVars;
depthVecSim = resSim.simData.srcVars.xVec;
dxSim = resSim.simData.srcVars.meshRes;

for i=1:5
    phiSimRaw(i,:,:,:) = sqrt(resSim.simData.phiLight{i}(:,121-envSize:121+envSize,121-envSize:121+envSize));
end
% Align MCX phi to Simulation Grid
phiSim = vAOSim.interpAlignReplPhi(phiSimRaw, depthVecSim, true, false, true);

%% Simulate Conv-based Fluence
resConvSP  = vAOSim.reconSP(abs(phiSim), true);
resConvHad = vAOSim.reconHad(abs(phiSim), true);

%% Load Measurements:
ml = measLoader();

uVarsML = ml.createUserVars();

uVarsML.projectPath  = "D:/MuEff/Uniform/Focused";
uVarsML.measNameTemp = "Phantom%d-Focused";
uVarsML.numMeas      = 5;

ml.setVars(uVarsML)
ml.loadMeas();
%
lowIdx  = 1;
highIdx = 235;

ml.cutPhi(lowIdx, highIdx);
%
c = 1400;
ml.fixSpeedOfSound(c)

dX = vAOSim.vars.dX;
ml.intAndAlign(dX);

%
tailIdx = 1370; % Post interpolation
ml.normToTail(tailIdx);

%
measData = ml.data;

%% 
res = vAOSim.matchMeasAndSpeckleSim(measData.phiRawAlignNorm,...
                              measData.xVecAligned,...
                              [],...
                              [],...
                              resConvSP.phiEnvReconNorm',...
                              phiSim(:,:,1,1));

%%





%% Simulate Speckle-based Fluence
clc; close all;
resSpeckle = vAOSim.speckleSim(abs(phiSim), false);
vAOSim.displaySpeckleRes(resSpeckle);

% Interpolate & Align Reconstructions:
speckSPAlign = vAOSim.interpAlignReplPhi(resSpeckle.phiReconSPNorm, resSpeckle.vars.xUS, false, true, false);
speckHadAlign = vAOSim.interpAlignReplPhi(resSpeckle.phiReconHadNorm, resSpeckle.vars.xUS, false, true, false);
%% 
vAOSim.matchMeasAndSpeckleSim(measData.phiRawAlignNorm,...
                              measData.xVecAligned,...
                              resSpeckle.phiReconHadNorm,...
                              resSpeckle.vars.xUS,...
                              [],...
                              phiSim(:,:,1,1));