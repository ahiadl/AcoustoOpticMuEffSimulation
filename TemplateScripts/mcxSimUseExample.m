close all
clear all
clc;

ms = mcxSim();

uVars = ms.createUserVars();

uVars.xLims = [-40, 40];
uVars.yLims = [-30, 30];
uVars.zLims = [0, 70];

uVars.res = 0.5;

uVars.objParams{1}.shape = 'rect';
uVars.objParams{1}.bb    = [-25, -25, 5;
                            -20,  -20, 10];
uVars.objParams{1}.op.mua = 30;
uVars.objParams{1}.op.mus = 4;
uVars.objParams{1}.op.g   = 0.55;
uVars.objParams{1}.op.ref = 1.34;

uVars.objParams{2}.shape  = 'sphere';
uVars.objParams{2}.center = [-15, 20, 15];
uVars.objParams{2}.radius = 5; 
uVars.objParams{2}.op.mua = 20;
uVars.objParams{2}.op.mus = 4;
uVars.objParams{2}.op.g   = 0.55;
uVars.objParams{2}.op.ref = 1.34;

uVars.objParams{3}.shape = 'cylinder';
uVars.objParams{3}.dir     = 'Z';
uVars.objParams{3}.center  = [-5,0];
uVars.objParams{3}.radius  = 5;
uVars.objParams{3}.longLim = [10,60];
uVars.objParams{3}.op.mua = 30;
uVars.objParams{3}.op.mus = 4;
uVars.objParams{3}.op.g   = 0.55;
uVars.objParams{3}.op.ref = 1.34;

uVars.objParams{4}.shape = 'layer';
uVars.objParams{4}.dir = 'X';
uVars.objParams{4}.depth = 0;
uVars.objParams{4}.op.mua = 1;
uVars.objParams{4}.op.mus = 4;
uVars.objParams{4}.op.g   = 0.55;
uVars.objParams{4}.op.ref = 1.34;

uVars.srcType = 'pencil';
uVars.srcPos  = [-40, 0, 35];
uVars.srcDir  = [1, 0, 0];
uVars.pattern = [];

uVars.planePoint1   = [];   %%
uVars.planePoint2   = [];   %%
uVars.patternWidth  = [];   %%
uVars.patternHeight = [];   %%
uVars.waistSize     = [];   %%
uVars.halfAngle     = [];   %%
uVars.lineEndPoint  = [];   %%
uVars.slitEndPoint  = [];   %%

uVars.tstart = 0;
uVars.tstep  = 1e-10;
uVars.tend   = 5e-9;

uVars.cw = true;

%Bkg optical properties:
uVars.mua = 0.004;  % 1/mm
uVars.mus = 0.4; % 1/mm
uVars.g   = 0.55;  % [-1,1]
uVars.ref = 1.34;  % 1<

uVars.nphoton = 1e7;

ms.setVars(uVars);

phi = ms.calcPhi();

xVec = ms.grid.xVec;
yVec = ms.grid.yVec;
zVec = ms.grid.zVec;

figure()
subplot(2,2,1)
hs=slice(yVec, xVec, zVec, log(phi), [-22.5], [-22.5], [7.5]);
set(hs,'linestyle','none');
axis tight equal
title("Rect"); xlabel("Y"); ylabel("X"); zlabel("Z");
subplot(2,2,2)
hs=slice(yVec, xVec, zVec, log(phi), [20],[-15],[15]);
set(hs,'linestyle','none');
axis tight equal
title("Sphere"); xlabel("Y"); ylabel("X"); zlabel("Z");
subplot(2,2,3)
hs=slice(yVec, xVec, zVec, log(phi), [0],[-5],[35]);
set(hs,'linestyle','none');
axis tight equal
title("Cylinder"); xlabel("Y"); ylabel("X"); zlabel("Z");
subplot(2,2,4)
hs=slice(yVec, xVec, zVec, log(phi), [0],[0],[35]);
set(hs,'linestyle','none');
axis tight equal
title("Layer"); xlabel("Y"); ylabel("X"); zlabel("Z");

figure()
subplot(2,2,1)
hs=slice(yVec, xVec, zVec, phi, [-22.5], [-22.5], [7.5]);
set(hs,'linestyle','none');
axis tight equal
title("Rect"); xlabel("Y"); ylabel("X"); zlabel("Z");
subplot(2,2,2)
hs=slice(yVec, xVec, zVec, phi, [20],[-15],[15]);
set(hs,'linestyle','none');
axis tight equal
title("Sphere"); xlabel("Y"); ylabel("X"); zlabel("Z");
subplot(2,2,3)
hs=slice(yVec, xVec, zVec, phi, [0],[-5],[35]);
set(hs,'linestyle','none');
axis tight equal
title("Cylinder"); xlabel("Y"); ylabel("X"); zlabel("Z");
subplot(2,2,4)
hs=slice(yVec, xVec, zVec, phi, [0],[0],[35]);
set(hs,'linestyle','none');
axis tight equal
title("Layer"); xlabel("Y"); ylabel("X"); zlabel("Z");


%% plane
ms = mcxSim();

uVars = ms.createUserVars();

uVars.xLims = [-30, 30];
uVars.yLims = [-30, 30];
uVars.zLims = [0, 40];

uVars.res = 1;

uVars.tstart = 0;
uVars.tstep  = 1e-10;
uVars.tend   = 5e-9;

uVars.cw = true;

uVars.numOfObjects = 1; %TODO: maybe should be in geometry?

uVars.mua = 0.04;   % 1/mm
uVars.mus = 26.66; % 1/mm
uVars.g   = 0.8;   % [0,1]
uVars.ref = 1.34;  % 1<

uVars.nphoton = 1e7;

uVars.srcType = 'planar';
uVars.srcPos  = [-30, 0, 20];
uVars.srcDir  = [1, 0, 0];

uVars.plane = 'XY';
uVars.plane = 'XZ';
uVars.plane = 'YZ';

uVars.planeSize = [20,20];

ms.setVars(uVars);
tic
phi = ms.calcPhi();
toc

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

%% Pattern
clc
close all
clear all

% Load Pattern:
ia = illuminationAnalysis();
[pattern, vars]=ia.analyse("31-Oct-2022 12-39-50-Illumination-NoPDMS.mat");

uVars.pattern = pattern.mask;
% close all

% pattern=[...
% 0 0 0 0 0 0 0 0 0 0 0
% 0 1 1 0 0 0 0 0 1 1 0
% 0 0 0 1 1 0 1 1 0 0 0
% 0 0 0 0 0 1 0 0 0 0 0
% 0 0 0 1 1 0 1 1 0 0 0
% 0 1 1 0 0 0 0 0 1 1 0
% 0 0 0 0 0 0 0 0 0 0 0
% 0 0 1 0 0 0 0 0 1 0 0
% 0 1 0 0 0 0 0 0 0 1 0
% 0 1 0 0 0 0 0 0 0 1 0
% 0 1 0 0 0 0 0 0 0 1 0
% 0 0 1 1 1 1 1 1 1 0 0
% 0 0 0 0 0 0 0 0 0 0 0
% 0 1 1 1 1 1 1 1 1 1 0
% 0 0 0 1 0 0 0 0 0 0 0
% 0 0 0 0 1 1 0 0 0 0 0
% 0 0 0 1 0 0 0 0 0 0 0
% 0 1 1 1 1 1 1 1 1 1 0
% 0 0 0 0 0 0 0 0 0 0 0];

%%
clc
% close all
ms = mcxSim();

uVars = ms.createUserVars();

uVars.xLims = [0, 70];
uVars.yLims = [-30, 30];
uVars.zLims = [-40, 40];

uVars.res = 0.5;

uVars.objParams{1}.shape = 'rect';
uVars.objParams{1}.bb    = [  5,  -5, -5;
                             10,  5, 5];
uVars.objParams{1}.op.mua = 100;
uVars.objParams{1}.op.mus = 4;
uVars.objParams{1}.op.g   = 0.55;
uVars.objParams{1}.op.ref = 1.34;

uVars.mua = 0.0046;   % 1/mm
uVars.mus = 4; % 1/mm
uVars.g   = 0.55;   % [0,1]
uVars.ref = 1.34;  % 1<

uVars.tstart = 0;
uVars.tstep  = 1e-10;
uVars.tend   = 5e-9;

uVars.cw = true;

uVars.nphoton = 1e7;
uVars.srcType = 'pattern';
uVars.pattern = pattern;
uVars.pattern = pattern.mask;
uVars.srcPos  = [0, 0, 0];
uVars.srcDir  = [0, 1, 0];

% uVars.srcPos  = [0, 0, 0];
% uVars.srcDir = [0, 0, 1];

% uVars.patternVec1 = 0:0.5:0.5*((size(uVars.pattern,1)-1));
% uVars.patternVec2 = 0:0.5:0.5*((size(uVars.pattern,2)-1));
% uVars.patternVec1 = 0:1:1*((size(uVars.pattern,1)-1));
% uVars.patternVec2 = 0:1:1*((size(uVars.pattern,2)-1));

uVars.patternVec1 = vars.scanVec - min(vars.scanVec);
uVars.patternVec2 = vars.disc1Vec - min(vars.disc1Vec);

% uVars.plane = 'XZ';
uVars.plane = 'YZ';
% uVars.plane = 'XY';

% uVars.plane = 'custom';
% uVars.planePoint1   = [19, 0, 0];
% uVars.planePoint2   = [0, 11, 0];

uVars.patternSizeFactor = 1;

ms.setVars(uVars);
%%
phi = ms.calcPhi();
phi(isnan(phi)) = 1e-10;

%%
xVec = ms.grid.xVec;
yVec = ms.grid.yVec;
zVec = ms.grid.zVec;

sizeG = ms.grid.size;

figure()
subplot(2,2,1)
imagesc(yVec, xVec, squeeze(phi(:,:,1)))
xlabel("Y[mm]"); ylabel("X[mm]");
axis tight equal
subplot(2,2,2)
imagesc(zVec, xVec, log(squeeze(phi(30,:,:))))
xlabel("Z[mm]"); ylabel("Y[mm]");
axis tight equal
subplot(2,2,3)
hs=slice(yVec, xVec, zVec, log(phi), [yVec(1), yVec(end)], [xVec(1), xVec(end)], [zVec(1), zVec(end)]);
set(hs,'linestyle','none');
axis tight equal
xlabel("Y"); ylabel("X"); zlabel("Z");
subplot(2,2,4)
hs=slice(yVec, xVec, zVec, phi, [0],[0],[0]);
set(hs,'linestyle','none');
axis tight equal
xlabel("Y"); ylabel("X"); zlabel("Z");

% figure()
% hs=slice(yVec, xVec, zVec, log(phi), [0],[0],[0]);
% set(hs,'linestyle','none');
% axis tight equal
% title("Layer"); xlabel("Y"); ylabel("X"); zlabel("Z");
