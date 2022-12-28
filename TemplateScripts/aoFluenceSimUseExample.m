close all
clear all
clc

aoSim = aoFluenceSim();

%%
clc
close all
uVarsAOSim = aoFluenceSim.createUserVars();

uVarsAOSim.usPath  = "Transducer Pressure Field\AO_Transducer.mat";
uVarsAOSim.srcPath = "31-Oct-2022 12-39-50-Illumination-NoPDMS.mat";
uVarsAOSim.detPath = "31-Oct-2022 16-56-07-Reflection-NoPDMS.mat";

uVarsAOSim.phiPath = [];
uVarsAOSim.loadPhi = false;

uVarsAOSim.numOfPhantom = 5;
uVarsAOSim.geometry = 'Measured';
uVars.srcParam      = [];

uVarsAOSim.simulator = 'mcx'; % 'mcx', 'toast'

uVarsAOSim.muaVec = [0.046, 0.092, 0.184, 0.368, 0.736]/10;
uVarsAOSim.mus    = [0.9];   %
uVarsAOSim.g      = [0.555]; %
uVarsAOSim.ref    = [1.34];  %

uVarsAOSim.nphoton = 1e8; %1e8

uVarsAOSim.meshSize  = [60,60,60]; % lenX, lenY, lenZ
uVarsAOSim.meshRes   = 0.25; % [mm]

uVarsAOSim.usAx = 'X';
uVarsAOSim.usDir = -1; % [1/-1]
uVarsAOSim.usFocalDist = 30; %distance of the US focus from illumination plane
uVarsAOSim.usDataType = '3D';

uVarsAOSim.savePath = '.';
uVarsAOSim.saveFlag = false;
uVarsAOSim.simName  = 'FocusedTrans';

aoSim.setVars(uVarsAOSim);
aoSim.simAndAnalyse();


%%
% res = load("13-Nov-2022 14-22-49-FocusedTrans.mat");
% 
% phiSrc   = normMatf(res.simData.phiSrc{1});
% phiDet   = normMatf(res.simData.phiDet{1});
% phiLight = normMatf(res.simData.phiLight{1});
% 
% xVec = res.simData.srcVars.xVec;
% yVec = res.simData.srcVars.yVec;
% zVec = res.simData.srcVars.zVec;
% 
% figure()
% subplot(3,3,1)
% imagesc(yVec, xVec, log(squeeze(phiSrc(:,:,121))))
% axis tight equal
% xlabel("Y[mm]")
% ylabel("X[mm]")
% subplot(3,3,4)
% imagesc(zVec, yVec, log(squeeze(phiSrc(121,:,:))))
% axis tight equal
% xlabel("Z[mm]")
% ylabel("Y[mm]")
% subplot(3,3,7)
% hs=slice(log(phiSrc), [1,241],[1, 241],[1 241]);
% set(hs,'linestyle','none');
% axis tight equal
% xlabel("Y[mm]")
% ylabel("X[mm]")
% zlabel("Z[mm]")
% 
% subplot(3,3,2)
% imagesc(yVec, xVec, log(squeeze(phiDet(:,:,121))))
% axis tight equal
% xlabel("Y[mm]")
% ylabel("X[mm]")
% subplot(3,3,5)
% imagesc(zVec, yVec, log(squeeze(phiDet(121,:,:))))
% axis tight equal
% xlabel("Z[mm]")
% ylabel("Y[mm]")
% subplot(3,3,8)
% hs=slice(log(phiDet), [1,241],[1, 241],[1 241]);
% set(hs,'linestyle','none');
% axis tight equal
% xlabel("Y[mm]")
% ylabel("X[mm]")
% zlabel("Z[mm]")
% 
% subplot(3,3,3)
% imagesc(yVec, xVec, log(squeeze(phiLight(:,:,121))))
% axis tight equal
% xlabel("Y[mm]")
% ylabel("X[mm]")
% subplot(3,3,6)
% imagesc(zVec, yVec, log(squeeze(phiLight(121,:,:))))
% axis tight equal
% xlabel("Z[mm]")
% ylabel("Y[mm]")
% subplot(3,3,9)
% hs=slice(log(phiLight), [1,241],[1, 241],[1 241]);
% set(hs,'linestyle','none');
% axis tight equal
% xlabel("Y[mm]")
% ylabel("X[mm]")
% zlabel("Z[mm]")
