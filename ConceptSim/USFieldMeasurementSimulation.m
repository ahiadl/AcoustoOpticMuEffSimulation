close all
clear all
clc

c = 1500;

% Define Signal:
f=1e6;
dts = 10e-9;
tVecSig=0:dts:(1/f-dts);
sig = sin(2*pi*f*tVecSig);

% Hydropphone Translation span):
spaceLen = length(sig) *15;
dx = c*dts;
xVec = (0:1:(spaceLen-1))*dx;

% Sampling rate and Sampled Signal Length:
Ts = 20e-6;
tsVec = 0:dts:(Ts-dts);
depthVec = tsVec*c;

% Pulse Geometrical Quantities:
pulseIdx = length(sig);
pulseLen = dx*pulseIdx;

% Create Envelope Function:
envelopeField = 22.5e3-xVec*1e6;

% Scan Geometry:
scanSpan = pulseIdx * 10;
sMin = pulseLen;
sMinIdx = pulseIdx; % in xVec
scanVec = (0:1:(scanSpan-1))*dx+sMin;

%Scan temporal properties:
delayT   = sMin/c;
delayIdx = delayT/dts;

% Sampled Signal model
sigS = zeros(1, length(tsVec));
sigS(1:pulseIdx) = sig;

figure();
subplot(1,2,1)
plot(xVec*1e3, envelopeField)
subplot(1,2,2)
plot(tsVec*1e6, sigS)

%% Simulate scan:
measMat = zeros(length(scanVec), length(tsVec));

for i = 1:length(scanVec)    
    measMat(i, :) = circshift(sigS, pulseIdx+i)* envelopeField(i+pulseIdx);
end

% figure();
% plot(circshift(sigS, pulseIdx+i))

figure();
subplot(1,2,1)
imagesc(depthVec*1e3, scanVec*1e3, measMat)
axis tight equal
xlabel('distanceFromDetector')
ylabel('detector Position')
subplot(1,2,2)
plot(tsVec, measMat(1,:)); hold on
plot(tsVec, measMat(250,:))
plot(tsVec, measMat(500,:))
plot(tsVec, measMat(750,:))
plot(tsVec, measMat(1000,:))


%% Reconstruct Space:
idx = 500;
x0 = xVec(idx);

spaceMat = zeros(length(scanVec), length(tsVec)+length(scanVec)-1);

for i=1:length(scanVec)
   spaceMat(i,i:i+length(tsVec)-1) = measMat(i,:);
end

%%
close all
figure()
subplot(2,2,1)
imagesc(spaceMat)
subplot(2,2,2)
imagesc(measMat)
subplot(2,2,3)
plot(spaceMat(:,600)); hold on
plot(measMat(:,400))