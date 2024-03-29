 function [] = analyzeSensitivityDistribution(dataparms)
    %% Last Updated 2/24/2021 JM
    %Todo: Implement a filter for bad cells in the initial data processing
        %Idea for filter: remove responses below -2 and above 2
    %Todo: Use keita's method for bootstrapping confidence interval on CDF
    %points. Bootstrap from individual responses instead of single cell average

    %This script is meant to be used to analyze data output by
    %"FretDataExtractionAndProcessing" by K.K. 

    %System parameters can be changed in 'parameters.m'.
    %Path to datasets can be selected in 'datasets.m'. Plotting parameters must
    %also be selected in this document. 

    %% Add necessary functions to the PATH
    close all; 
%     clearvars -except dataparms;
    addpath('./Functions')
    addpath('./Data')
    addpath('./Outputs')
    parms = parameters(); %Struct containing system parameters
%     showPlot = parms.showPlot;

    %% Datasets: Load 
    % dataparms = datasets(); %Load data parameters. Not needed anymore, since
    % dataparms is provided upstream
%     backConc = dataparms.backConc;
%     xlabels = dataparms.xlabels;
%     Lplot = dataparms.Lplot;
    Files = dataparms.Files;
    OutFolder = dataparms.OutFolder; 
    concLevels = dataparms.concLevels;

    %% Make output destination
    OutputDest = ['./Outputs/', OutFolder];
    mkdir(OutputDest)
    OutputDest

    %% Initial data processing/reformatting
    EfretData = combineFretFiles(Files, concLevels); %Combine files into one dataset
    n = length(concLevels); %Number of stimulus levels
    m = parms.m; %Number of responses per stim level

    %Reformat data into a struct with n matrices with m measurements each
    doseData = reorganizeData(EfretData, n, m, OutputDest);

    %Calculate response-amplitude for each individual stimulus
    doseData = calcResponseAmp2(doseData, n, parms);
    badCells = find([doseData.a0] < 0.15);
    doseData(badCells) = [];
   

    %Calculate population-average and single-cell dose-response data
    [popAvgData, SCAvgData] = calcPopAvgDoseResponse2(m, doseData, parms, dataparms);

    
    
    %% For QC, get x and y position of each cell
%     %After looking at all the data, I am convinced that there is no spatial
% %     dependence on the responses. Good news.
    
%     AvgResponseByCell = squeeze(median(popAvgData.allRData, 1));
%     xPos = [];
%     yPos = [];
%     for ii = 1:length(EfretData)
%         xPos = [xPos, EfretData(ii).xPos];
%         yPos = [yPos, EfretData(ii).yPos];
%     end
%     
%     respno = 5;
%     figure();
%     subplot(2, 1, 1)
%     plot(xPos, AvgResponseByCell(respno,:), 'o')
%     xlabel("x")
%     subplot(2, 1, 2)
%     plot(yPos, AvgResponseByCell(respno,:), 'o')
%     xlabel("y")

    %% Distribution of dynamic ranges
    %Find EfretMax - EfretMin and (EfretMax - EfretMin) / EfretMin
    FretRangeData = findRangeDistribution(EfretData, OutputDest);

    %% Plot FRET Time-series
%     close all
%     plotFRETTimeSeries(EfretData, dataparms, parms, OutputDest);
       
    %% Plot Distribution of Responses for Each Level Over Time
    close all
    plotResponseDistribution(EfretData, dataparms, parms, OutputDest, m)

    %% Fit hill function to population-average dose-response data
    %Collect dose-response data, and optimal parameters for hill function
    HillPlotData = fitHillFunction(popAvgData, SCAvgData, dataparms, parms, OutputDest);

    %% Fit Hill Function to Single-cell dose-response data
%     fitSCHillFunction(popAvgData, SCAvgData, dataparms, parms, OutputDest);
    
    %% Fit CDF
    close all
    CDFPlotData = fitLogNormCDF(SCAvgData, popAvgData, dataparms, parms, OutputDest);

    %% Save Information for Recreating Plots Elsewhere
    save([OutputDest, 'plotData.mat'], "HillPlotData", "CDFPlotData", "FretRangeData");

    %% End
    % close all; 
end
