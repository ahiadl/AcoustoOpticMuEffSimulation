close all
clear all
clc;

%%
usa = usAnalysis();

%% 3D
uVars = usa.createUserVars();

uVars.usDataPath = 'AO_Transducer.mat';
uVars.usDataType = '3D';
uVars.intFactor = [1,1,1];
usa.setVars(uVars);
res = usa.analyse();
%%
close all
% usa.cutFocalData();
usa.interpolateFocalEnvelopes();
% usa.createSpatialData();
% usa.calcPulseWidth();
% usa.cutPulses();
% close all;
% usa.displayResults();

%% slice

uVars.usDataPath = '..\Measurements\Transducer Pressure Field\02-May-2022 16-36-41-PA-Flat-XAxis-FullScan-1.25MHz.mat';
uVars.usDataType = '2D';
usa.setVars(uVars);
res = usa.analyse();


%% tmp
xScan = load("..\Measurements\Transducer Pressure Field\02-May-2022 16-36-41-PA-Flat-XAxis-FullScan-1.25MHz.mat");
xVec = xScan.csVars.scanVec;
tVec = xScan.csVars.tVec;
dx = xVec(2) - xVec(1);
dt = tVec(2) - tVec(1);
data = xScan.resCs;

% P2P
p2p = max(data, [], 5) - min(data, [], 5);

figure();
subplot(1,2,1)
plot(xVec, p2p);
subplot(1,2,2)
plot(squeeze(data(1,:)));
% hold on;
% plot(squeeze(data(100,:)));
% plot(squeeze(data(500,:)));

% Focus idx
[~, focusIdx] = max(p2p);

% Speed of sound
sigF = squeeze(data);
for i=1:20
    [~, idxT(i)] = max(abs(sigF(focusIdx +(i-1),:)));
end
dTpeak =  abs(mean(tVec(idxT(2:end)) - tVec(idxT(1:end-1))));
c = dx*1e-3/dTpeak;

%Calc Depth Vec
depthVec = c * tVec * 1e3; %[mm];
dDepth = abs(depthVec(2) - depthVec(1)); 

% Calc Delay in focus
idx = focusIdx;
idxTrig = 405;
sigMin  = data(idx, :);
sigG    = gradient(sigMin);
startTIdxOfSig = find(sigG(idxTrig+1:end) > 0.01*max(sigG),1) + idxTrig;
startTIdxOfSig = 6577;
figure();
plot(squeeze(data(focusIdx,:)));
delayAtIdx = tVec(startTIdxOfSig);
FocalLen = c*delayAtIdx*1e3;










