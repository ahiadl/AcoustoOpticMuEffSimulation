close all
clear all
clc;
%%
usl = usLoader();

uVars = usLoader.createUserVars();
uVars.paths{1} = ".\AcoustoOpticSystem\Measurements\Transducer Pressure Field\Calibrated\Focused\1.25MHz-Depth-Axis.mat";
uVars.paths{2} = ".\AcoustoOpticSystem\Measurements\Transducer Pressure Field\Calibrated\Focused\1.25MHz-Transversal-Axis.mat";

uVars.intFactor      = [2.5,2.5];
uVars.usAtZero       = false;
uVars.useExtSpeed    = false;
uVars.extC           = 1487;
uVars.pulseSizeIdx   = 120;
uVars.envPeakSize    = 10;
uVars.displayResults = false;

resFocused = usl.loadAndAnalyse(uVars);

usl.displayResults();

dirname = "./AcoustoOpticSystem/Measurements/Transducer Pressure Field/Calibrated";
filename = sprintf("%s/Analysis - FocusedAOITransducer-1.25MHz.mat", dirname);
% save(filename, '-Struct', 'resFocused', '-v7.3');
%%
close all

usl = usLoader();

uVars = usLoader.createUserVars();
uVars.paths{1} = ".\AcoustoOpticSystem\Measurements\Transducer Pressure Field\Calibrated\UnFocused\1.25MHz-Depth-Axis.mat";
uVars.paths{2} = ".\AcoustoOpticSystem\Measurements\Transducer Pressure Field\Calibrated\UnFocused\1.25MHz-Transversal-Axis.mat";

uVars.intFactor      = [2.5,2.5];
uVars.usAtZero       = false;
uVars.useExtSpeed    = false;
uVars.extC           = 1487;
uVars.pulseSizeIdx   = 120;
uVars.envPeakSize    = 15;
uVars.displayResults = false;

resUnFocused = usl.loadAndAnalyse(uVars);

usl.displayResults();

dirname = "./AcoustoOpticSystem/Measurements/Transducer Pressure Field/Calibrated";
filename = sprintf("%s/Analysis - UnFocusedAOITransducer-1.25MHz.mat", dirname);
% save(filename, '-Struct', 'resUnFocused', '-v7.3');

%% 
close all
idxFocused = resFocused.grid.int.axialFocalIdx;
idxUnFocused = resUnFocused.grid.int.axialFocalIdx;
focalPulseFocused = resFocused.data.int.depth(idxFocused, :);
focalPulseUnFocused = resUnFocused.data.int.depth(idxUnFocused, :);
axFocused = resFocused.grid.int.depthVec;
axUnFocused = resUnFocused.grid.int.depthVec;

pulseFNorm = focalPulseFocused / max(focalPulseFocused);
pulseUFNorm = focalPulseUnFocused / max(focalPulseUnFocused);
pulseUFNorm = circshift(pulseUFNorm, -122);

figure();
plot(axFocused, pulseFNorm); hold on
plot(axUnFocused, pulseUFNorm); 
xlabel("Depth [mm]")
ylabel("Acoustic Intensity")
xlim([88, 100])


legend("Focused", "Flat")