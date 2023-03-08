close all;
clear all;
clc;

%%
ia = illuminationAnalysis();
ia.analyse("31-Oct-2022 12-39-50-Illumination-NoPDMS.mat");
ia.analyse("31-Oct-2022 16-56-07-Reflection-NoPDMS.mat");

%%
close all

ia = illuminationAnalysis();
ia.analyse("31-Jan-2023 08-59-25-Transmission-Uniform-NoSpacer.mat");
ia.analyse("31-Jan-2023 11-57-49-Reflection-Uniform-NoSpacer.mat");
ia.analyse("31-Jan-2023 11-14-00-Transmission-Point-NoSpacer.mat");
ia.analyse("31-Jan-2023 10-56-30-Reflection-Point-NoSpacer.mat");

