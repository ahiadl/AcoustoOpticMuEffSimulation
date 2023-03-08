close all
clear ll
clc;

%%
close all
clc

ml = measLoader();

uVarsML = ml.createUserVars();

uVarsML.projectPath  = 'D:/MuEff/Uniform/GradedScatteringSet/UnFocused';
uVarsML.measNameTemp = "Phantom-%d";
uVarsML.numMeas      = 5;

uVarsML.idxLow        = 1;
uVarsML.idxHigh       = 236;
uVarsML.phiHighToLow  = true;
uVarsML.noiseIdxs     = 1:30;
uVarsML.loadNew       = true;

ml.setVars(uVarsML)
data = ml.loadMeas();
%
ml.resetCurData();
c = []; % use [] to disable SOS correction
ml.fixSpeedOfSound(c);

ml.interp(0.01)
ml.align();

ml.displayResults();

%%
fields{1,1} = 'extPeakIdx';
fields{1,2} = 95;
fields{2,1} = 'extNoiseIdx';
fields{2,2} = 30:60;
fields{3,1} = 'extSNRDef';
fields{3,2} = true;
fields{4,1} = 'noiseStd';
fields{4,2} = [];
pivot= 95;
env = 5;
offset = pivot - env-1;

for i=1:5
    curPhi = data.resMeas(i).phiNorm;
    [~, I] = max(curPhi(pivot-env : pivot + env));
    fields{1,2} = I + offset;
    resNew(i) = ml.reCalcMeasAndSave(fields, i);
end

figure();
for i=1:5
    plot(resNew(i).phiNorm); hold on
end
