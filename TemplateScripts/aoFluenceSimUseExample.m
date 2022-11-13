close all
clear all
clc

aoSim = aoFluenceSim();

%%
clc
close all
uVarsAOSim = aoFluenceSim.createUserVars();

uVarsAOSim.usPath  = "..\Measurements\Transducer Pressure Field\AO_Transducer.mat";
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

uVarsAOSim.nphoton = 1e7; %1e8

uVarsAOSim.meshSize  = [60,60,60]; % lenX, lenY, lenZ
uVarsAOSim.meshRes   = 0.25; % [mm]

uVarsAOSim.usAx = 'X';
uVarsAOSim.usDir = -1; % [1/-1]
uVarsAOSim.usFocalDist = 30; %distance of the US focus from illumination plane
uVarsAOSim.usDataType = '3D';

aoSim.setVars(uVarsAOSim);
aoSim.simAndAnalyse();


phiLight = aoSim.data.phiLight;
%%
% aoSim.data.phiLight = zeros(size(aoSim.data.phiLight));
% aoSim.data.phiLight(1,121,121) = 1;
aoSim.data.phiLight = phiLight;
aoSim.data.phiAO = aoSim.cnv.calcConv(aoSim.data.phiLight, aoSim.data.srcVars);

% figure();
% imagesc(squeeze(aoSim.data.phiLight(:,121,:)))

figure();
subplot(2,2,1)
imagesc(log(normMatf(squeeze(aoSim.data.phiAO(:,129, :)))))
subplot(2,2,2)
imagesc(squeeze(aoSim.data.phiAO(:,129, :)))
subplot(2,2,3)
plot(squeeze(aoSim.data.phiAO(129,129, :)))