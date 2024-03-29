%% Explanation
% This is meant to generate supplemental figures for figure 2. The main
% supplement is a plot of the K1/2 CVs with 95 confidence intervals, so
% it's possible to tell if distributions are similar or not.

close all; clear;
addpath('./Functions')

% MeAsp dose responses
Files1 = ["./Outputs/MeAspDR_OneConcentrationSets/", "./Outputs/100MeAspBackMeAspDR/",...
    "./Outputs/1SerBackMeAspDR/", "./Outputs/100MeAsp1SerBackMeAspDR/"];

Files2 = ["./Outputs/0BackSerDose/", "./Outputs/1SerBackSerDose/",...
    "./Outputs/100MeAspBackSerDose/", "./Outputs/100MeAsp1SerBackSerDR/"];

%% Plot statistics of CV curves in figure 2
figure()

[CV1, means1, stds1, CVerr1, CVlowerbound1, CVupperbound1, meanerr1, ...
    meanlowbound1, meanupperbound1] = extract_cv_from_files(Files1);
subplot(2, 2, 1)
background = ["None", "100\muM meAsp", "1\muM Ser", "1\muM Ser, 100\muM meAsp"];
x = categorical(background);
CVPlot = scatter(x, CV1, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
hold on
CVError = errorbar(x, CV1, CV1-CVlowerbound1, CVupperbound1-CV1, "LineStyle", "none", "color", "k"); 
% CVError = errorbar(x, CV, CVerr, CVerr, "LineStyle", "none", "color", "k"); 
ylabel("meAsp K_{1/2} CV")

subplot(2, 2, 3)
MeanPlot = scatter(x, means1, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
hold on
% MeanError = errorbar(x, means, meanerr, meanerr, "LineStyle", "none", "color", "k"); 
MeanError = errorbar(x, means1, means1-meanlowbound1, meanupperbound1-means1, "LineStyle", "none", "color", "k");
set(gca, 'yscale', 'log')
ylabel("Mean meAsp K_{1/2}")


[CV2, means2, stds2, CVerr2, CVlowerbound2, CVupperbound2, meanerr2,...
    meanlowbound2, meanupperbound2] = extract_cv_from_files(Files2);
subplot(2, 2, 2)
background = ["None", "1\muM Ser", "100\muM meAsp", "100\muM meAsp 1\muM Ser"];
x = categorical(background);
CVPlot = scatter(x, CV2, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
hold on
CVError = errorbar(x, CV2, CV2-CVlowerbound2, CVupperbound2-CV2, "LineStyle", "none", "color", "k"); 
% CVError = errorbar(x, CV, CVerr, CVerr, "LineStyle", "none", "color", "k"); 
ylabel("Ser K_{1/2} CV")

subplot(2, 2, 4)
MeanPlot = scatter(x, means2, "MarkerEdgeColor", "k", "MarkerFaceColor", "k");
hold on
% MeanError = errorbar(x, means, meanerr, meanerr, "LineStyle", "none", "color", "k"); 
MeanError = errorbar(x, means2, means2-meanlowbound2, meanupperbound2-means2, "LineStyle", "none", "color", "k");
set(gca, 'yscale', 'log')
ylabel("Mean Ser K_{1/2}")

saveas(gcf, './Outputs/Figure2_supplement.svg')

%% Plot comparison of meAsp K1/2 distributions in the presence of maltose

Files1 = {["./Outputs/MeAspDR_OneConcentrationSets/", ...
    "./Outputs/1MaltosebackMeAspDose/"]};

% Color Palette To Color Each plot. Rows organized the same way
colrs = {["#4E4E4E", "#6aa84f"]};

lnstyles = {["-", "-"]};

CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));


%Position of the plots organized in the same order as Files1
plotInds = [1, 3];

xlims = [[10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3]];

xlabels = ["[meAsp] \muM", "", "[meAsp] \muM", ""];

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