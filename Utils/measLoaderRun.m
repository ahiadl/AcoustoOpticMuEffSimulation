close all
clear all
clc

ml = measLoader();

uVars = ml.createUserVars();

uVars.projectPath  = "D:/MuEff/Uniform/Focused";
uVars.measNameTemp = "Phantom%d-Focused";
uVars.numMeas      = 5;


ml.setVars(uVars)
ml.loadMeas();
%%
lowIdx  = 1;
highIdx = 235;

ml.cutPhi(lowIdx, highIdx);
%%
dX = 0.1492;
ml.intAndAlign(dX);

%%
tailIdx = 1370; % Post interpolation
ml.normToTail(tailIdx);

%%
data = ml.data;