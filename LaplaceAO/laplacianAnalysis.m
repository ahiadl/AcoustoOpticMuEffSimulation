close all
clear all
clc;

ms = mcxSim();

uVars = ms.createUserVars();

uVars.xLims = [0, 70];
uVars.yLims = [-30, 30];
uVars.zLims = [-40, 40];

uVars.res = 0.5;

uVars.objParams{1}.shape  = 'sphere';
uVars.objParams{1}.center = [50, 0, 0];
uVars.objParams{1}.radius = 10; 
uVars.objParams{1}.op.mua = 0.05;
uVars.objParams{1}.op.mus = 0.8;
uVars.objParams{1}.op.g   = 0.55;
uVars.objParams{1}.op.ref = 1.37;

uVars.srcType = 'pencil';
uVars.srcPos  = [1, 1, 1];
uVars.srcDir  = [1, 0, 0];
uVars.pattern = [];

uVars.tstart = 0;
uVars.tstep  = 1e-10;
uVars.tend   = 5e-9;

uVars.cw = true;

%Bkg optical properties:
uVars.mua = 0.004;  % 1/mm
uVars.mus = 0.8; % 1/mm
uVars.g   = 0.55;  % [-1,1]
uVars.ref = 1.37;  % 1<

uVars.nphoton = 1e10;

ms.setVars(uVars);
vol = ms.cfg.vol;
phi = ms.calcPhi();

muEffBkg = sqrt(3*0.004*(0.004+0.8*(1-0.55)));
muEffAbs = sqrt(3*0.050*(0.050+0.8*(1-0.55)));

xVec = ms.grid.xVec;
yVec = ms.grid.yVec;
zVec = ms.grid.zVec;

figure()
subplot(1,2,1)
hs=slice(yVec, xVec, zVec, log(phi), [0], [5], [0]);
set(hs,'linestyle','none');
axis tight equal
title("Rect"); xlabel("Y"); ylabel("X"); zlabel("Z");
subplot(1,2,2)
hs=slice(yVec, xVec, zVec, phi, [0], [5], [0]);
set(hs,'linestyle','none');
axis tight equal
title("Rect"); xlabel("Y"); ylabel("X"); zlabel("Z");

%%
[Y, X, Z] = meshgrid(yVec, xVec, zVec);

%Filtering:
%------------
%--- Gaussian
phiFilt = smooth3(phi,'gaussian') ;
%--- Median:
% phiFilt = medfilt3(phi) ;
%--- No Filter:
% phiFilt = phi;

% First Derivative:
%------------------
[gY, gX, gZ] = gradient(phiFilt);
grad = sqrt(gY.^2 + gX.^2 +gZ.^2);
grad2 = grad./phiFilt;

% Second Derivative:
%-------------------
gX = gX/uVars.res; gY = gY/uVars.res; gZ=gZ/uVars.res;
lap = divergence(Y, X, Z, gY, gX, gZ);

% Full Recon:
% -----------
recon = sqrt(abs(lap./phiFilt));

reconFilt = smooth3(recon,'gaussian', [7,7,7]) ;
absPhi = reconFilt(52:69, 61, 69:91);
reconMuEff = mean(absPhi(:))

figure()
subplot(2,2,1); imagesc(zVec, xVec, squeeze(vol(:,61,:)))
title("Grad"); xlabel("Z[mm]"); ylabel("X[mm]");
colorbar; axis tight equal
subplot(2,2,2); imagesc(log(squeeze(phi(:,61,:))))
title("Phi"); xlabel("Z[mm]"); ylabel("X[mm]");
colorbar; axis tight equal
subplot(2,2,3); imagesc(log(squeeze(grad2(:,61,:))))
title("Grad"); xlabel("Z[mm]"); ylabel("X[mm]");
colorbar; axis tight equal
subplot(2,2,4); imagesc(squeeze(reconFilt(:, 61,:)));
title("Laplacian"); xlabel("Z[mm]"); ylabel("X[mm]");
colorbar; axis tight equal

figure();
subplot(1,2,1)
imagesc(zVec, xVec, squeeze(reconFilt(:, 61,:)));
title("Laplacian"); xlabel("Z[mm]"); ylabel("X[mm]");
colorbar; axis tight equal
subplot(1,2,2)
plot(xVec, squeeze(recon(:, 61,81)));
title("Laplacian"); xlabel("Z[mm]"); ylabel("\mu_{eff}");

figure();
subplot(1,2,1)
imagesc(squeeze(reconFilt(:, 61,:)));
title("Laplacian"); xlabel("Z[mm]"); ylabel("X[mm]");
colorbar; axis tight equal
subplot(1,2,2)
imagesc(squeeze(reconFilt(71, :,:)));
title("Laplacian"); xlabel("Z[mm]"); ylabel("Y[mm]");
colorbar; axis tight equal

 
%% 