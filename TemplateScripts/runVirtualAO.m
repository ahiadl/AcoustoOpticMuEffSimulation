close all;
clear all;
clc;

vAOSim = virtualAOSim();
vAOSim.loadUSPulse();
%%
close all

uVars = vAOSim.createUserVars();

uVars.debug = false;
uVars.displayDebug = false;
uVars.debugTime = false;

uVars.N   = 251;
uVars.fUS = 1.25e6;
uVars.pulseType = 'measured'; %delta, measured, 

uVars.spacerLen = 6.4; %mm
uVars.spacerMaterial = 'Water';

uVars.usFocalPointDistFromInterface = 30; %[mm] 

uVars.muEffVec      = [0.074; 0.106; 0.152; 0.219; 0.323];

vAOSim.setVars(uVars)
simVars = vAOSim.getVars();

vAOSim.extractUSVars();
vAOSim.useCustomUSParams(true)
vAOSim.calcSimDimension();
vAOSim.createSpace();
vAOSim.createMathematicalFluence();
vAOSim.alignAndInterpUS();
vAOSim.createAcousticInterface();
vAOSim.createPulses();
vAOSim.buildSPMatrix();
vAOSim.buildHadMatrix();
vAOSim.buildHadInvMat();
% vAOSim.reconAll(vAOSim.phiMath);
% vAOSim.displayResults();
% resSpeckle = vAOSim.speckleSim(vAOSim.phiMath, true);
%% 
% clc
% close all
% resSpeckle = vAOSim.speckleSim(vAOSim.phiMath, true);

%% Load MCX Simulations:
resSim    = load("13-Nov-2022 14-22-49-FocusedTrans.mat");

phiSimRaw = zeros(5,241);
varsSim = resSim.simVars;
depthVecSim = resSim.simData.srcVars.xVec;
dxSim = resSim.simData.srcVars.meshRes;

for i=1:5
    phiSimRaw(i,:) = resSim.simData.phiLight{i}(:,121,121);
end
%% Align MCX phi to Simulation Grid
phiSim = vAOSim.alignExternalPhi(phiSimRaw, depthVecSim, true);

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
c = vAOSim.vars.us.c;
ml.fixSpeedOfSound(c)

dX = vAOSim.vars.dX;
ml.intAndAlign(dX);

%
tailIdx = 1370; % Post interpolation
ml.normToTail(tailIdx);

%
measData = ml.data;

%% Simulate Speckle-based Fluence
clc; close all;
resSpeckle = vAOSim.speckleSim(abs(phiSim), false);
vAOSim.displaySpeckleRes(resSpeckle);
%% 
vAOSim.matchMeasAndSpeckleSim(measData.phiRawAlignNorm,...
                              measData.xVecAligned,...
                              resSpeckle.phiReconSPNorm,...
                              resSpeckle.vars.xUS,...
                              resConvSP.phiEnvReconNorm',...
                              phiSim);