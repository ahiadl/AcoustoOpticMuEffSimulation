% close all
% clear all
clc

pd = PhantomDesigner();

uVars.wl   = 785; %[nm]
uVars.Vph  = 200;
uVars.musP = 0.4*ones(1,5); %[mm^-1]
uVars.mua  = 0.00258 * (2.^[1:5]); %[mm^-1]
uVars.numLayers = 1;
uVars.numPh = 5;
uVars.inkType = 'RoyalTalens';
uVars.layerThick = 56; %[mm]
%%
% clc
pd.createPhantomSet(uVars);

%% Graded Absorption:
clc;
uVars.musP = 0.8; %[mm^-1]
uVars.mua  = 2.6e-3 * 2.^(1:5); %[mm^-1]
uVars.mua  = 2.6e-3 * 2 * (1+(0:4)*0.5); %[mm^-1]
pd.gradedAbsorption(uVars);

%% Graded Scattering:1200
clc;
uVars.musP = 0.02*2.^(1:5); %[mm^-1]
uVars.mua  = 2*2.58e-3;
pd.gradedScattering(uVars);

%%
D = 54.8;
L = 80;

vol = (D^2*L)*1e-3;
mo100w = 83.8;
moWeight = vol/100*mo100w;


% muaNigBase = 

