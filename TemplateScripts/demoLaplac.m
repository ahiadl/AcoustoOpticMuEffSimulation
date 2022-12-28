%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCXLAB - Monte Carlo eXtreme for MATLAB/Octave by Qianqina Fang
%
% In this example, we demonstrate how to use sub-pixel resolution 
% to represent the problem domain. The domain is consisted of a 
% 6x6x6 cm box with a 2cm diameter sphere embedded at the center.
%
% This file is part of Monte Carlo eXtreme (MCX) URL:http://mcx.sf.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cfg;

% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E'); 

% format: [mua(1/mm) mus(1/mm) g n]
cfg.prop=[0     0   1   1    ;   % medium 0: the environment
          0.004 0.8 0.5 1.37 ;   % medium 1: cube
          0.050 0.8 0.5 1.37];   % medium 2: spherical inclusion

muEffBkg = sqrt(3*0.004*(0.004+0.8*(1-0.5)));
muEffAbs = sqrt(3*0.050*(0.050+0.8*(1-0.5)));

% time-domain simulation parameters
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-10;

% GPU thread configuration
cfg.autopilot=1;
cfg.gpuid=1;

cfg.isreflect=1; % enable reflection at exterior boundary
cfg.isrefint=1;  % enable reflection at interior boundary too

% define the source position
cfg.srcpos=[60,60,0]+1;
cfg.srcdir=[0 0 1];

% define a 1cm radius sphere within a 6x6x6 cm box with a 0.5mm resolution
dim=120;
[xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);
dist=(xi-40).^2+(yi-40).^2+(zi-60).^2;
cfg.vol=ones(size(xi));
cfg.vol(dist<400)=2;
cfg.vol=uint8(cfg.vol);

cfg.unitinmm = 0.5; % define the pixel size in terms of mm
cfg.nphoton = 1e9;  % you need to simulate 8x photons to get the same noise

[f2,det2]=mcxlab(cfg);

%%
load('pointPhi.mat');
phi = res.phi;
X = res.X;
Y = res.Y;
Z = res.Z;

%%
phi = sum(f2.data,4)*cfg.tstep;
X = (xi-40)*0.5; Y = (yi-40)*0.5; Z = zi *0.5;
%%
% phiFilt = medfilt3(phi,[3 3 3]);
phiFilt = smooth3(phi,'gaussian') ;
phiFilt = phi;
[gY, gX, gZ] = gradient(phiFilt);
grad = sqrt(gY.^2 + gX.^2 +gZ.^2);
grad2 = grad./phiFilt;
gX = gX/0.5; gY = gY/0.5; gZ=gZ/0.5;
lap = divergence(X, Y, Z, gX, gY, gZ);

recon = sqrt(abs(lap./phiFilt));

figure()
subplot(2,2,1)
imagesc(squeeze(cfg.vol(:,40,:)))
title("Grad")
colorbar
axis tight equal
subplot(2,2,3)
imagesc(log(squeeze(grad2(:,40,:))))
title("Grad")
colorbar
axis tight equal
subplot(2,2,2)
imagesc(log(squeeze(phi(:,40,:))))
title("Phi")
colorbar
axis tight equal
subplot(2,2,4)
imagesc(squeeze(recon(:, 40,:)));
title("laplacian")
colorbar
axis tight equal

figure();
plot(squeeze(recon(40,40,:)))

res.phi = phi;
res.X = X;
res.Y = Y;
res.Z = Z;
res.op = cfg.prop;
res.muEffBkg = muEffBkg;
res.muEffAbs = muEffAbs;
res.vol = cfg.vol;
save('pointPhi.mat', 'res', '-v7.3');