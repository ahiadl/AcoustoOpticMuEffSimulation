close all
clear all
clc;

ms = mcxSim();

uVars = ms.createUserVars();

uVars.xLims = [-30, 30];
uVars.yLims = [-30, 30];
uVars.zLims = [0, 40];

uVars.res = 1;

uVars.srcType = 'pencil';
uVars.srcPos  = [30, 30, 1];
uVars.srcDir  = [0, 0, 1];
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

uVars.numOfObjects = 1; %TODO: maybe should be in geometry?

uVars.mua = 0.04;   % 1/mm
uVars.mus = 26.66; % 1/mm
uVars.g   = 0.8;   % [-1,1]
uVars.ref = 1.34;  % 1<

uVars.nphoton = 1e7;

ms.setVars(uVars);
tic
phi = ms.calcPhi();
toc

figure()
imagesc(log(squeeze(phi(30,:,:))))

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
close all
ms = mcxSim();

uVars = ms.createUserVars();

uVars.xLims = [0, 60];
uVars.yLims = [-30, 30];
uVars.zLims = [-30, 30];

uVars.res = 1;

uVars.mua = 0.0046;   % 1/mm
uVars.mus = 26.66; % 1/mm
uVars.g   = 0.8;   % [0,1]
uVars.ref = 1.34;  % 1<

uVars.tstart = 0;
uVars.tstep  = 1e-10;
uVars.tend   = 5e-9;

uVars.cw = true;

uVars.numOfObjects = 1; %TODO: maybe should be in geometry?

uVars.nphoton = 1e7;
uVars.srcType = 'pattern';
uVars.pattern = pattern;
uVars.pattern = pattern.mask;
uVars.srcPos  = [0, 0, 0];
uVars.srcDir  = [1, 0, 0];

% uVars.srcPos  = [0, 0, 0];
% uVars.srcDir = [0, 0, 1];

% uVars.patternVec1 = 0:0.5:0.5*((size(uVars.pattern,1)-1));
% uVars.patternVec2 = 0:0.5:0.5*((size(uVars.pattern,2)-1));
% uVars.patternVec1 = 0:1:1*((size(uVars.pattern,1)-1));
% uVars.patternVec2 = 0:1:1*((size(uVars.pattern,2)-1));

uVars.patternVec1 = vars.scanVec-min(vars.scanVec);
uVars.patternVec2 = vars.disc1Vec - min(vars.disc1Vec);

% uVars.plane = 'XZ';
uVars.plane = 'YZ';
% uVars.plane = 'XY';

% uVars.plane = 'custom';
% uVars.planePoint1   = [19, 0, 0];
% uVars.planePoint2   = [0, 11, 0];

uVars.patternSizeFactor = 1;

ms.setVars(uVars);

tic
phi = ms.calcPhi();
toc

phi(isnan(phi)) = 1e-10;

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



