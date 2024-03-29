%% Control script for running sensitivity distribution analysis
%This version of the code is designed to be more friendly when analyzing
%multiple datasets at the same time. The analyses are kept independent
%until it comes time to plot dose-response and cumulative distribution
%functions

clear all;
close all;
% dataparms = datasets_ShortSat();
% dataparms = datasets();
dataparms = datasets_ForPublication();

for ii = [2]%length(dataparms)
%     dataparms = dataparms(i);
    analyzeSensitivityDistribution(dataparms(ii));
end

%% Get nCells from a dataset
dataparms = datasets_ForPublication();
addpath('./Functions')
addpath('./Data')
addpath('./Outputs')
ii = 45;
Files = dataparms(ii).Files;
concLevels = dataparms(ii).concLevels;
EfretData = combineFretFiles(Files, concLevels); %Combine files into one dataset
disp(length(EfretData))


