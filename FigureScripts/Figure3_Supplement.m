%% Explanation
%This script is to generate figure S3 of the diversity tuning paper. The
%goal for this figure is to display statistical analysis of 2 datasets. On
%the top will be analysis of the meAsp Khalf distribution in the presence
%of various competitors. On the bottom will be analysis of all the
%Tar-binding ligands Khalf distributions in the absence of any background. 

%TODO: Make ticks bigger

close all; clear;
addpath('./Functions')

%% Load data and declare plot parameters for meAsp khalf with competitors

%Files are organized by row:
    %Row 1: top left
    %Row 2: bottom left
    %Row 3: top right
    %Row 4: bottom right
Files1 = {["./Outputs/MeAspDR_OneConcentrationSets/",  "./Outputs/100GluBackMeAspDR/", "./Outputs/1mMGluBackMeAspDose/",...
    "./Outputs/1LAsp_meAspDR/",...
    "./Outputs/10LAspBackMeAspDose/",...
    "./Outputs/1MaltosebackMeAspDose/"]};

% Color Palette To Color Each plot. Rows organized the same way
colrs = {["#4E4E4E", "#008080", "#008080", "#e02459",  "#e02459", "#FF934F"]};

lnstyles = {["-", "--", "-", "--", "-", "-"]};

CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));


%%=============================================
%% Plot K1/2 distributions all on the same axes
%%=============================================
%Establish grid: Should be 4x2

%Position of the plots organized in the same order as Files1
plotInds = [1, 3];

xlims = [[10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3]];

xlabels = ["[L] \muM", "", "[meAsp] \muM", ""];

figure()
for ii = 1:length(Files1)
    FilesTemp = Files1{ii};
    colr = colrs{ii};
    xlimits = xlims(ii, :);
    lnstyl = lnstyles{ii};

    for jj = 1:length(FilesTemp)
        %Load Data
        load([convertStringsToChars(FilesTemp(jj)), 'plotData.mat'])
        Lplot = CDFPlotData.Lplot;

        %Generate sample curves for error bars
        rnd = CDFPlotData.MCMCSamples; %Dimensions: [mu, sigma, A (not used)]
        

        %Plot PDF on top
        subplot(4, 2, plotInds(ii, 1))
        hold on
        LplotPDF = log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot));
        y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
        plot(LplotPDF, y ./ max(y), 'color', colr(jj), 'Linewidth', 1.8,...
            'LineStyle', lnstyl(jj))
%         ylabel(['normalized', newline, 'pdf'])
        xlim([log(xlimits(1)), log(xlimits(2))])
        set(gca, 'XTick', [])
        if plotInds(ii, 2) /2 == round(plotInds(ii, 2)/2)
            set(gca, 'YTick', [])
        end
        if plotInds(ii, 2) /2 ~= round(plotInds(ii, 2)/2)
            ylabel('pdf')
        end
        if(jj == length(FilesTemp))
            set(gca, 'Position', get(gca, 'Position') + [0, -0.02, 0, -0.1])
        end
        

        %Plot CDF on bottom
        curveSamples = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples(j, :) = CDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr = prctile(curveSamples, 97.5);
        curveNegErr = prctile(curveSamples, 2.5);

        subplot(4, 2, plotInds(ii, 2))
        hold on
        shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(jj))
        errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr,...
            CDFPlotData.CDFPointPosErr, 'o', 'color', colr(jj));
        C(jj) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), ...
            'color', colr(jj), 'Linewidth', 1.8, 'LineStyle', lnstyl(jj));
        xlabel(xlabels(ii))
        xlim(xlimits)
        
        set(gca, 'xscale', 'log', 'XTick', [10^-2, 10^-1, 10^0, 10^1, 10^2])
        if plotInds(ii, 2) /2 ~= round(plotInds(ii, 2)/2)
            ylabel('P(K_{1/2} < [L])')
        else
            set(gca, 'YTick', [])
        end
        
    end
end

set(gcf, 'Position', [384,54,844,898])

%%======================================
%% Plot mean, std, and CV for all curves
%%======================================
CVFiles = Files1{1};
[CV, means, stds, CVerr, CVlowerbound, CVupperbound, meanerr, meanlowbound, meanupperbound] = extract_cv_from_files(CVFiles);
% Order of files: [0, 100Glu, 1mMglu, 1LAsp, 10LAsp]

subplot(4, 2, 2)
background = ["None", "100\muM Glu", "1mM Glu", "1\muM LAsp", "10\muMLAsp", "1\muM Maltose"];
x = categorical(background);
CVPlot = scatter(x, CV, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
hold on
CVError = errorbar(x, CV, CV-CVlowerbound, CVupperbound-CV, "LineStyle", "none", "color", "k"); 
% CVError = errorbar(x, CV, CVerr, CVerr, "LineStyle", "none", "color", "k"); 
ylabel("meAsp K_{1/2} CV")
ylim([0, 1.7])

subplot(4, 2, 4)
MeanPlot = scatter(x, means, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
hold on
% MeanError = errorbar(x, means, meanerr, meanerr, "LineStyle", "none", "color", "k"); 
MeanError = errorbar(x, means, means-meanlowbound, meanupperbound-means, "LineStyle", "none", "color", "k");
set(gca, 'yscale', 'log', 'YTick', [1, 2, 4, 10, 20])
ylabel("Mean meAsp K_{1/2}")

saveas(gcf, './Outputs/Figure3_supplement.svg')




%% ==========================================
%% Khalf distribution for Tar-binding ligands
%% ==========================================
%% Load data and declare plot parameters for tar-ligand khalf distributions

%Files are organized by row:
    %Row 1: top left
    %Row 2: bottom left
    %Row 3: top right
    %Row 4: bottom right
Files1 = {["./Outputs/MeAspDR_OneConcentrationSets/",  "./Outputs/0BackLAspDose/", ...
    "./Outputs/0BackGluDose/",...
    "./Outputs/0backMaltoseDose/",...
    "./Outputs/1MaltosebackMaltoseDose/"]};

% Color Palette To Color Each plot. Rows organized the same way
colrs = {["#4E4E4E", "#e02459", "#008080",  "#FF934F", "#FF934F"]};

lnstyles = {["-", "-", "-", "-", "--"]};

CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));


%%=============================================
%% Plot K1/2 distributions all on the same axes
%%=============================================
%Establish grid: Should be 4x2

%Position of the plots organized in the same order as Files1
plotInds = [1, 3];

xlims = [[10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3]];

xlabels = ["[L] \muM", "", "[meAsp] \muM", ""];

figure()
for ii = 1:length(Files1)
    FilesTemp = Files1{ii};
    colr = colrs{ii};
    xlimits = xlims(ii, :);
    lnstyl = lnstyles{ii};

    for jj = 1:length(FilesTemp)
        %Load Data
        load([convertStringsToChars(FilesTemp(jj)), 'plotData.mat'])
        Lplot = CDFPlotData.Lplot;

        %Generate sample curves for error bars
        rnd = CDFPlotData.MCMCSamples; %Dimensions: [mu, sigma, A (not used)]
        

        %Plot PDF on top
        subplot(4, 2, plotInds(ii, 1))
        hold on
        LplotPDF = log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot));
        y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
        plot(LplotPDF, y ./ max(y), 'color', colr(jj), 'Linewidth', 1.8,...
            'LineStyle', lnstyl(jj))
%         ylabel(['normalized', newline, 'pdf'])
        xlim([log(xlimits(1)), log(xlimits(2))])
        set(gca, 'XTick', [])
        if plotInds(ii, 2) /2 == round(plotInds(ii, 2)/2)
            set(gca, 'YTick', [])
        end
        if plotInds(ii, 2) /2 ~= round(plotInds(ii, 2)/2)
            ylabel('pdf')
        end
        if(jj == length(FilesTemp))
            set(gca, 'Position', get(gca, 'Position') + [0, -0.02, 0, -0.1])
        end
        

        %Plot CDF on bottom
        curveSamples = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples(j, :) = CDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr = prctile(curveSamples, 97.5);
        curveNegErr = prctile(curveSamples, 2.5);

        subplot(4, 2, plotInds(ii, 2))
        hold on
        shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(jj))
        errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr,...
            CDFPlotData.CDFPointPosErr, 'o', 'color', colr(jj));
        C(jj) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), ...
            'color', colr(jj), 'Linewidth', 1.8, 'LineStyle', lnstyl(jj));
        xlabel(xlabels(ii))
        xlim(xlimits)
        
        set(gca, 'xscale', 'log', 'XTick', [10^-2, 10^-1, 10^0, 10^1, 10^2])
        if plotInds(ii, 2) /2 ~= round(plotInds(ii, 2)/2)
            ylabel('P(K_{1/2} < [L])')
        else
            set(gca, 'YTick', [])
        end
        
    end
end

set(gcf, 'Position', [384,54,844,898])

%%======================================
%% Plot mean, std, and CV for all curves
%%======================================
CVFiles = Files1{1};
[CV, means, stds, CVerr, CVlowerbound, CVupperbound, meanerr, meanlowbound, meanupperbound] = extract_cv_from_files(CVFiles);
% Order of files: [0, 100Glu, 1mMglu, 1LAsp, 10LAsp]

subplot(4, 2, 2)
background = ["meAsp", "L-Asp", "Glu", "Maltose", "1uM Maltose"];
x = categorical(background);
CVPlot = scatter(x, CV, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
hold on
CVError = errorbar(x, CV, CV-CVlowerbound, CVupperbound-CV, "LineStyle", "none", "color", "k"); 
% CVError = errorbar(x, CV, CVerr, CVerr, "LineStyle", "none", "color", "k"); 
ylabel("K_{1/2} CV")
ylim([0, 1.7])

subplot(4, 2, 4)
MeanPlot = scatter(x, means, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
hold on
% MeanError = errorbar(x, means, meanerr, meanerr, "LineStyle", "none", "color", "k"); 
MeanError = errorbar(x, means, means-meanlowbound, meanupperbound-means, "LineStyle", "none", "color", "k");
set(gca, 'yscale', 'log', 'YTick', [0.1, 1, 10, 100])
ylabel("Mean K_{1/2}")

saveas(gcf, './Outputs/Figure3_supplement_2.svg')

