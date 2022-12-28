%% Init Matlab
close all;
clear all;
clc;

%% Load US pulse:
if false
    usa = usAnalysis();

    uVars.usDataPath = 'AO_Transducer.mat';
    uVars.usDataType = '3D';
    uVars.intFactor = [1,1,1];
    usa.setVars(uVars);
    res = usa.analyse();
    usVars = usa.getVars();
    % usa.cutPulses()
    % usa.calcPulsesProfile();
    % usa.extractFocalSignal();
    
    resSave.focalProfile      = res.focalProfile;
    resSave.focalPulsesEnvCut = res.focalPulsesEnvCut;
    resSave.focalSig          = res.focalSig;
    resSave.focalSigEnv       = res.focalSigEnv;
    resSave.focalPulseRawAx   = usVars.depthVec;

    save("C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\Simulations\AcoustoOpticMuEffSimulation\analysedFocusedUS.mat", 'resSave', 'usVars', '-v7.3');
    
end
load("C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\Simulations\AcoustoOpticMuEffSimulation\analysedFocusedUS.mat");

focalProfileMat   = resSave.focalProfile;
focalPulsesEnvCut = resSave.focalPulsesEnvCut;
focalSig          = resSave.focalSig;
focalSigEnv       = resSave.focalSigEnv;
focalPulseRawAx   = resSave.focalPulseRawAx; 

%% US Extraction:
focalProfile     = squeeze(focalProfileMat(9,9,:));
focalProfileNorm = focalProfile/max(focalProfile);
depthProfile     = envelope(focalProfile, 25, 'peaks');
[maxVal, profileFocalPointIdx] = max(depthProfile);
depthProfileNorm = depthProfile/max(depthProfile);

pulse = squeeze(focalPulsesEnvCut(9,9, usVars.focalIndIntAx,:))';
pulseNorm = normMatf(pulse);
pulseAx = usVars.pulseVec;

focalPulseSig = flip(focalSig(649:649+149-1));
focalPulseSigNorm = focalPulseSig/max(focalPulseSig);

profileDepthVecRaw = usVars.pulseStartPosVec;
usFocalLen = profileDepthVecRaw(profileFocalPointIdx);

figure();
subplot(1,2,1)
plot(pulseAx, pulse);
title("Pulse Envelope")
xlabel("X[mm]")
ylabel("Normalized Pressure")
subplot(1,2,2)
plot(profileDepthVecRaw, focalProfileNorm);hold on
plot(profileDepthVecRaw, depthProfileNorm)
title("US Beam Axial Profile");
xlabel("X[mm]")
ylabel("Normalized Pressure")

%% Simulation Parameters:
debug = false;
c = round(usVars.c);
N =  251;

fUS = 1.25e6;

dX = usVars.dAx;
singleCycleLen = (1/1.25e6)*c*1e3;
singleCycleIdxRaw = singleCycleLen/dX;
singleCycleIdx = round(singleCycleIdxRaw);

reconSize = N*singleCycleIdx;

usFocalPointDistFromInterface = 30; %[mm]
spacerLen = 6.4*2;

numOfPhantoms = 5;
muEffVec = [0.074; 0.106; 0.152; 0.219; 0.323];

%% Create Mathematical Fluence:
xRaw = (0:1:reconSize-1)*dX;
x = xRaw - (usFocalLen+usFocalPointDistFromInterface); % mm

spacerX   = 0:dX:spacerLen +dX;
spacer    = zeros(5,length(spacerX));

phiLenIdx = length(x);

phiMath1 = exp(-muEffVec.*abs(x));
phiMath2 = exp(-muEffVec.*abs(x-spacerX(end)-dX));

phiMath = zeros(length(muEffVec), length(x));
phiMath(:, x<=0)                   = phiMath1(:, x<=0);
phiMath(:, x>0 & x<= spacerX(end)) = 0;
phiMath(:, x>spacerX(end))         = phiMath2(:, x>spacerX(end));

figure();
subplot(1,2,1)
plot(x, phiMath')
title("Mathematical Fluence Profile - with Reflection")
xlabel("X[mm]")
ylabel("Fluence [AU]")
subplot(1,2,2)
plot(x, log(phiMath'))
title("Mathematical Fluence Profile - with Reflection")
xlabel("X[mm]")
ylabel("Fluence [Log]")

%% US Alignment & Adaptations
[~, focalPointIdx] = min(abs(x+usFocalPointDistFromInterface));
usFocalPos = x(focalPointIdx);

profileDepthVec = profileDepthVecRaw - usFocalLen + usFocalPos;
profileDepthVecArt = [x(1), profileDepthVec, x(end)];
focalProfileInt = normMatf(interp1(profileDepthVecArt, [0, depthProfile', 0], x, 'pchip'));

% figure();
% subplot(1,2,1)
% plot(pulseAx, pulse);
% subplot(1,2,2)
% plot(profileDepthVec, focalProfileNorm);hold on
% plot(x, focalProfileInt)

%% Test Pulse
% pulseInt = zeros(1, pulseIntLenIdx);
% pulseInt (end-11:end) = sin((0:11)*pi/6);

%% Convolution Example:

% kernel = [zeros(1,9), 1];
% kernel = linspace(0,1,10);
% f = linspace(1,100, 100);
% f = [zeros(1,49), 1, zeros(1,50)];
% convRes = conv(f,flip(kernel));
% figure(); 
% subplot(2,2,1); plot(flip(kernel));
% subplot(2,2,2); plot(f)
% subplot(2,1,2); plot(convRes);

%% Acoustic Reflection Profile:
zPDMS = 1.048e6;  %[Ns * m^-3]
zWater = 1.494e6; %[Ns * m^-3]
Gamma = abs((zPDMS - zWater) / (zPDMS + zWater));
usRefProfile(x<=0) = 1;
usRefProfile(x>0)  = Gamma;

figure()
hold on
plot(x, phiMath);
plot(x, focalProfileInt)
plot(x, usRefProfile);
xlabel("X[mm]")
ylabel("AU")

%% Naive Simulation:
measPhiNaive = phiMath.*focalProfileInt.*usRefProfile;

figure();
for i=1:5
    ax(i) = subplot(2,3,i);
    hold on
    plot(log(phiMath(i,:)))
    plot(real(log(measPhiNaive(i,:))))
    xlabel("X[mm]")
    ylabel("Fluence [AU]")
    title(sprintf("Math. Phantom: %d", i))
end
linkaxes(ax);

%% SP propagation matrix Considering Reflection Profile, US Profile and Pulse Width:
sVecInt    = zeros(1,reconSize);
sVecInt(1) = 1;

if debug
    
    pulseSP = [zeros(1,length(pulse)-1), 1];
    pulseSP = pulseNorm;
    pulseSigSP = focalPulseSigNorm;
    
    pulseSP    = pulseNorm(end-singleCycleIdx+1:end);
    pulseSigSP = focalPulseSigNorm(end-singleCycleIdx+1:end);
else
    pulseSP = pulse;
    pulseSigSP = focalPulseSigNorm;
end

padSize  = length(pulseSP)-1;
endIdx   = reconSize+padSize;
startIdx = padSize+1;

if debug
    figure()
    subplot(3,2,1); h1 = plot(zeros(1,length(sVecInt)));   title("CurSqnc"); ylim([-1,1]);
    subplot(3,2,2); h2 = plot(zeros(1,2*length(sVecInt))); title("curSqncSpatial2"); ylim([-1,1]);
    subplot(3,2,3); h3 = plot(zeros(1,2*length(sVecInt))); title("transmissionSP"); ylim([-1,1]);
    subplot(3,2,4); h4 = plot(zeros(1,length(sVecInt)));   title("transMatSP"); ylim([-1,1]);
    subplot(3,2,5); h5 = plot(zeros(1,2*length(sVecInt))); title("transmissionSP"); ylim([-1,1]);
    subplot(3,2,6); h6 = plot(zeros(1,length(sVecInt)));   title("transMatSP"); ylim([-1,1]);
end

transMatSP = zeros(reconSize);
transMatSigSP = zeros(reconSize);
for i = 1:reconSize
    curSqnc          = circshift(sVecInt, (i-1));
    curSqncSpatial   = curSqnc.*focalProfileInt.*usRefProfile;
    if debug
%         curSqncSpatial = curSqnc;
    end
    transmissionSP  = conv(curSqncSpatial, pulseSP, 'full');
    transMatSP(i,:) = transmissionSP(startIdx:endIdx);
    
    transmissionSigSP  = conv(curSqncSpatial, pulseSigSP, 'full');
    transMatSigSP(i,:) = transmissionSigSP(startIdx:endIdx);
    
    if debug
        set(h1, 'YData', curSqnc);
        set(h2, 'YData', curSqncSpatial);
        set(h3, 'YData', transmissionSP);
        set(h4, 'YData', transMatSP(i,:));
        set(h5, 'YData', transmissionSigSP);
        set(h6, 'YData', transMatSigSP(i,:));
        drawnow();
    end
end

transMatSPNorm = transMatSP/max(transMatSP(:));
transMatSPSigNorm  = transMatSigSP/max(transMatSigSP(:));

figure();
subplot(1,2,1)
imagesc(transMatSPNorm);
axis tight equal
colorbar
subplot(1,2,2)
imagesc(transMatSPSigNorm);
axis tight equal
colorbar

%% SP Env with Mathematical Phantom

% Envelope Reconstruction
phiMathSPRecon = transMatSPNorm * phiMath';
phiMathSPReconNorm = phiMathSPRecon ./ max(phiMathSPRecon,[],1);

figure();
for i=1:5
    ax(i) = subplot(2,3,i);
    hold on
    plot(log(phiMath(i,:)))
    plot(real(log(phiMathSPReconNorm(:,i))))
    xlabel("X[mm]")
    ylabel("Fluence [AU]")
    title(sprintf("Phantom: %d", i))
end
linkaxes(ax);

% Signal Reconstruction (meaningless)
phiMathSPRecon = transMatSPSigNorm * phiMath';
phiMathSPReconNorm = phiMathSPRecon ./ max(phiMathSPRecon,[],1);

figure();
for i=1:5
    ax(i) = subplot(2,3,i);
    hold on
    plot(log(phiMath(i,:)))
    plot(real(log(phiMathSPReconNorm(:,i))))
    xlabel("X[mm]")
    ylabel("Fluence [AU]")
    title(sprintf("Phantom: %d", i))
end
linkaxes(ax);

%% Build Hadamard Matrix with US Considerations:
close all

if debug
    NHad=11;
    singleCycleIdxHad = 50;
    reconSizeHad = NHad*singleCycleIdxHad;
%  pulseHad = pulse(end-singleCycleIdx+1:end);
%     pulseHad = [zeros(1,singleCycleIdx-1),1];
%     pulseHad = 11:1:70;
    pulseHad = pulseNorm(end-singleCycleIdx+1:end);
    pulseSigHad = focalPulseSigNorm(end-singleCycleIdx+1:end);
    
    
    
%     NHad = N;
%     singleCycleIdxHad = singleCycleIdx;
%     reconSizeHad = reconSize;
%     pulseHad = [zeros(1, length(pulse)-1), 1];
%     pulseHad = pulseNorm;
%     pulseSigHad = focalPulseSigNorm;
else
    NHad = N;
    singleCycleIdxHad = singleCycleIdx;
    reconSizeHad = reconSize;
    pulseHad = pulseNorm;
    pulseSigHad = focalPulseSigNorm;
end

sMat   = createSMatrix(NHad);
sVec   = sMat(1, :);

sVecInt          = zeros(singleCycleIdxHad,NHad);
sVecInt(end,:)   = sVec;
sVecInt          = flip(sVecInt(:)');

padSize = length(pulseHad)-1;

transFullLen = 3*reconSizeHad+padSize;

startIdx = reconSizeHad+padSize+1;
endIdx   = 2*reconSizeHad+padSize;

if debug
    figure()
    ax(1) = subplot(3,2,1); h1 = plot(zeros(1,length(sVecInt)));   title("CurSqnc");
    ax(2) = subplot(3,2,2); h2 = plot(zeros(1,2*length(sVecInt))); title("curSqncSpatial2");
    ax(3) = subplot(3,2,3); h3 = plot(zeros(1,2*length(sVecInt))); title("transmission");
    ax(4) = subplot(3,2,4); h4 = plot(zeros(1,length(sVecInt)));   title("transMat");
    ax(5) = subplot(3,2,5); h5 = plot(zeros(1,2*length(sVecInt))); title("transmission");
    ax(6) = subplot(3,2,6); h6 = plot(zeros(1,length(sVecInt)));   title("transMat");
    
    hold(ax(2), 'on');
    hold(ax(3), 'on');
    hold(ax(5), 'on');
    idxVec = 1:transFullLen;
    frame  = (idxVec >= startIdx) & (idxVec <= endIdx);
    h7 = plot(ax(2), frame);
    h8 = plot(ax(3), frame);
    h9 = plot(ax(5), frame);
end

transMat = zeros(reconSizeHad);

for i = 1:reconSizeHad
    curSqnc        = circshift(sVecInt, (i-1));
    if debug
%         curSqncSpatial = curSqnc;
        curSqncSpatial = curSqnc.*focalProfileInt.*usRefProfile;
    else
        curSqncSpatial = curSqnc.*focalProfileInt.*usRefProfile;
    end
    curSqncSpatial2 = [curSqncSpatial, curSqncSpatial, curSqncSpatial];
    
    transmission   = conv(curSqncSpatial2, pulseHad, 'full');
    transMat(i,:)  = transmission(startIdx:endIdx);
    
    transmissionSig   = conv(curSqncSpatial2, pulseSigHad, 'full');
    transSigMat(i,:)  = transmissionSig(startIdx:endIdx);
    
    if debug
        set(h1, 'YData', curSqnc);
        set(h2, 'YData', curSqncSpatial2);
        set(h3, 'YData', transmission);
        set(h4, 'YData', transMat(i,:));
        set(h5, 'YData', transmissionSig);
        set(h6, 'YData', transSigMat(i,:));
        drawnow();
    end
end

transMatHadNorm = transMat/ max(transMat(:));
transMatHadSigNorm = transSigMat/ max(transSigMat(:));

figure();
subplot(1,2,1)
imagesc(transMatHadNorm)
axis tight equal
title("Envelope Multiplexing Mat");
subplot(1,2,2)
imagesc(transMatHadSigNorm)
axis tight equal
title("Signal Multiplexing Mat");

%% Build clean Naive inverse matrix
sMatInv = inv(sMat);
sVecInv = sMatInv(1,:);

sVecInvSqnc = zeros(singleCycleIdx,N);
sVecInvSqnc(end,:) = sVecInv;
sVecInvSqnc = sVecInvSqnc(:)';

sMatInvSqnc = zeros(reconSize);
% 
% if debug
%     figure();
%     h = imagesc(sMatInvSqnc);
% end

for i=1:reconSize
    sMatInvSqnc(i,:) = circshift(sVecInvSqnc, i);
%     if debug
%         set(h, 'CData', sMatInvSqnc);
%         drawnow();
%     end
end
%% Apply Hadamard and Reconstruct Mathematical Fluence
close all;
phiClean          = phiMath';
% phiClean          = (1:550)';
phiReconClean     = sMatInvSqnc * (transMatHadNorm * phiClean);
phiReconCleanNorm = phiReconClean./max(phiReconClean,[],1);

figure();
subplot(1,2,1)
plot(x, phiReconCleanNorm)
xlabel("X[mm]")
ylabel("Fluence [AU]")
xlim([-60, 30])
title("Mathematical Fluence Profile")
subplot(1,2,2)
plot(x, log(abs(phiReconCleanNorm)))
xlabel("X[mm]")
ylabel("Fluence [log]")
xlim([-60, 30])
title("Mathematical Fluence Profile")

if debug
    figure();
    for i=1:5
        subplot(2,3,i)
        hold on
        plot(log(phiMath(i,:)))
        plot(log(abs(phiReconCleanNorm(:,i))))
        xlabel("X[mm]")
        ylabel("Fluence [AU]")
        title(sprintf("Phantom: %d", i))
    end
else  
    figure();
    for i=1:5
        subplot(2,3,i)
        hold on
        plot(x, log(phiMath(i,:)))
        plot(x, log(abs(phiReconCleanNorm(:,i))))
        xlabel("X[mm]")
        ylabel("Fluence [AU]")
        title(sprintf("Phantom: %d", i))
        xlim([-60,30]);
        ylim([-10,0]);
    end
end



%% Align Simulation Data
depthVecSimInt = depthVecSim(1) : dX: depthVecSim(end);
phiSim1SideInt = zeros(5,length(depthVecSimInt));
for i=1:5
    phiSim1SideInt(i,:) = interp1(depthVecSim, phiSimRaw(i,:), depthVecSimInt, 'pchip');
end

phiSim1SideIntFlip = flip(phiSim1SideInt,2);
phiSim2Side = [phiSim1SideIntFlip, spacer, phiSim1SideInt];
phiSim2SidePreInt = [zeros(5,1), phiSim2Side, zeros(5,1)];

depthVec2Side = (0:1:(size(phiSim2Side,2)-1))*dX;
xSim = depthVec2Side - depthVec2Side(length(phiSim1SideInt));
xSimPreInt = [x(1), xSim, x(end)];

phiSim = zeros(5, length(x));
for i=1:5
    phiSim(i,:) = interp1(xSimPreInt, phiSim2SidePreInt(i,:), x, 'pchip');
end

figure();
subplot(1,2,1)
plot(x, phiSim')
subplot(1,2,2)
plot(x, log(normMatf(phiSim,2)))

%% Apply Hadamard & Reconstruct
phiSimNorm = normMatf(phiSim, 2);

phiClean = phiSim';

phiReconClean = sMatInvSqnc * (transMatHadNorm * phiClean);
% phiReconCleanNorm = phiReconClean./max(phiReconClean, [], 1);
phiReconCleanNorm = normMatf(phiReconClean,1);
figure();
for i=1:5
    subplot(2,3,i)
    hold on
    plot(log(phiSimNorm(i,:)))
    plot(real(log(phiReconCleanNorm(:,i))))
    xlabel("X[mm]")
    ylabel("Fluence [AU]")
    title(sprintf("Phantom: %d", i))
end

figure();
hold on
plot(x, log(phiReconCleanNorm))
xlabel("X[mm]")
ylabel("Fluence [log]")
xlim([-60, 10])
ylim([-10,0])

figure();
for i=1:5
    subplot(2,3,i)
    hold on
    plot(x, log(phiSimNorm(i,:)))
    plot(x, log(phiReconCleanNorm(:,i)))
    xlabel("X[mm]")
    ylabel("Fluence [AU]")
    title(sprintf("Phantom: %d", i))
    xlim([-60,20]);
    ylim([-10,0]);
end

%% Speckle Sim:
close all
clear sig
Ephii = phiMath(1,:);

numIter = 8;

cl = 3e8;
n = 1.34;
lambda = 780e-9;
gamma = 0.001;

dt = usVars.tVec(2) -usVars.tVec(1); 
fs = 1/dt;

samplesPerPosPerSD = numIter*reconSize;
tVec = (0:1:numIter*reconSize-1)*dt;
fBarRaw = (0 : 1 : samplesPerPosPerSD-1) * (fs/samplesPerPosPerSD);

k0 = 2*pi/lambda; 
dnSP = gamma*transMatSPSigNorm;
dnHad = gamma*transSigMat;

sigSP  = zeros(reconSize, numIter);
sigHad = zeros(reconSize, numIter);
sigDM  = zeros(reconSize, numIter);

for i =1 :numIter
    phaseArg = repmat(2*pi*rand(1,reconSize), reconSize, 1);    
%     phaseArg = zeros(reconSize);

    wave  = repmat(exp(1i*n*k0.*x),reconSize, 1);
    phase = exp(1i*phaseArg);
    modUSSP = exp(1i*dnSP*k0.*x);
    modUSHad = exp(1i*dnHad*k0.*x);
    
    EiSP = sqrt(Ephii).*wave.*phase.*modUSSP;
    EiHad = sqrt(Ephii).*wave.*phase.*modUSHad;
    ItotSP = sum(EiSP,2).*conj(sum(EiSP,2));
    ItotHad = sum(EiHad,2).*conj(sum(EiHad,2));
    
    sigSP(:, i) = ItotSP;
    sigHad(:,i) = ItotHad;
    sigDM(:, i) = flip(sMatInvSqnc,2)*ItotHad;
end

A0 = reshape(sigSP',numIter,singleCycleIdx,N);
A1 = permute(A0, [2,1,3]);
A2 = reshape(A1, numIter*singleCycleIdx, 1, N);
A3 = permute(A2, [1,3,2]);

A4 = abs(fftshift(fft(A3,[],1))).^2;

figure();
subplot(1,2,1)
plot(ItotSP);
subplot(1,2,2)
plot(flip(sMatInvSqnc,2)*ItotHad)

%% Optimal Inversion Matrix:
close all

figure();
subplot(1,2,1)
hPs = plot(zeros(1,reconSize));
hTP = title("Plot number");
subplot(1,2,2)
hIm = imagesc(zeros(reconSize));

dmMat = zeros(reconSize,reconSize);

for i=1:reconSize
    if ~mod((i-1), 10)
    fprintf("Iter: %d/n", i);
    end
    curInvMat = circshift(sMatInvSqnc, (i-1),1);
    for j=1:reconSize
        curInvMat2 = circshift(curInvMat, (j-1),2);
        dmMat(j,:) = curInvMat2*ItotHad;
%         set(hPs, 'YData', squeeze(dmMat(i,j,:)))
%         set(hTP, 'String', num2str(j))
%         drawnow()
    end
    set(hIm, 'CData', squeeze(dmMat(:,:)))
    drawnow()
end



        


