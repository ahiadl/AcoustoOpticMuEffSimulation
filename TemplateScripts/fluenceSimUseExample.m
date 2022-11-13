close all
clear all
clc;

ia = illuminationAnalysis();
[pattern, vars]=ia.analyse("31-Oct-2022 12-39-50-Illumination-NoPDMS.mat");

fs = fluenceSim('mcx');

uVarsFS = fs.createUserVars();

uVarsFS.simulator = 'mcx'; % 'mcx', 'toast'
uVarsFS.mode      = 'Measured'; % twoFibers; Uniform; Measured
uVarsFS.meshMode  = 'slab'; % mcxOnly: 'slab', 'semiInf'
uVarsFS.meshSize  = [60,60,60]; % lenX, lenY, lenZ
uVarsFS.srcPos    = [0,0,0]; % [X,Y,Z] in mm; In case of two fibers this should be vector that specifices distance and angle
uVarsFS.srcSize   = []; % [mm]
uVarsFS.meshRes   = 1; % [mm]
uVarsFS.intFactor = 1; %
uVarsFS.srcDir    = [1,0,0]; % direction, as a normalized vector

uVarsFS.mua = [0.0001]; %
uVarsFS.mus = [26.66]; %
uVarsFS.g   = [0.8]; %
uVarsFS.ref = [1.34]; %

uVarsFS.loadPattern = false;

uVarsFS.pattern = pattern.mask;
uVarsFS.patternVec1 = vars.scanVec  - min(vars.scanVec);
uVarsFS.patternVec2 = vars.disc1Vec - min(vars.disc1Vec);

uVarsFS.nphoton = 1e7;

uVarsFS.plane = 'YZ';

fs.setVars(uVarsFS);
fs.config();
phi = fs.simulate();

figure()
subplot(2,2,1)
imagesc(log(squeeze(phi(:,:,1))))
imagesc(squeeze(phi(:,:,1)))
axis tight equal
subplot(2,2,2)
imagesc(log(squeeze(phi(30,:,:))))
axis tight equal
subplot(2,1,2)
hs=slice(phi, [1,61],[1, 61],[1 41]);
set(hs,'linestyle','none');
axis tight equal

