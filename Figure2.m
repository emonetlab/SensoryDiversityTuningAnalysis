%% Explanation
%This script is to generate figure 2 of the diversity tuning paper. The
%goal for this figure is to display 8 plots: in the top row will be the pdf
%and cdf of the K1/2 distribution for meAsp and serine in 0 background and
%a meAsp or serine background respectively. The bottom row will be the same
%but with the addition of either serine of meAsp, to show that the curves
%do not change. In inkscape, I will then add dotted lines connecting the
%top and bottom to show that the distributions are very similar. 

%TODO: Make ticks bigger

close all; clear;
addpath('./Functions')

%% Load data and declare plot parameters

%Files are organized by row:
    %Row 1: top left
    %Row 2: bottom left
    %Row 3: top right
    %Row 4: bottom right
Files1 = {["./Outputs/MeAspDR_OneConcentrationSets/", "./Outputs/100MeAspBackMeAspDR/"],...
        ["./Outputs/1SerBackMeAspDR/", "./Outputs/100MeAsp1SerBackMeAspDR/", "./Outputs/100MeAsp1SerBackMeAspDR/",...
        "./Outputs/MeAspDR_OneConcentrationSets/", "./Outputs/100MeAspBackMeAspDR/"],...
        ["./Outputs/0BackSerDose/", "./Outputs/1SerBackSerDose/"],...
        ["./Outputs/100MeAspBackSerDose/", "./Outputs/100MeAsp1SerBackSerDR/", "./Outputs/100MeAsp1SerBackSerDR/",...
        "./Outputs/0BackSerDose/", "./Outputs/1SerBackSerDose/"]};

% Color Palette To Color Each plot. Rows organized the same way
colrs = {["#4E4E4E", "#058ED9"],...
        ["#FF934F", "#FF934F", "#058ED9", "#4E4E4E", "#4E4E4E"],...
        ["#4E4E4E", "#FF934F"],...
        ["#058ED9", "#FF934F", "#058ED9", "#4E4E4E", "#4E4E4E"]};

lnstyles = {["-", "-"],...
        ["-",  "-", "--", ":", ":"],...
        ["-", "-"],...
        ["-", "-", "--", ":", ":"]};

showMarkers = {[1, 1],...
        [1,  1, 1, 0, 0],...
        [1, 1],...
        [1, 1, 1, 0, 0]};

lnwdths = {[1.8, 1.8],...
        [1.8,  1.8, 1.8, 1.2, 1.2],...s
        [1.8, 1.8],...
        [1.8, 1.8, 1.8, 1.2, 1.2]};

backgrColors = {[1 1 1],...
                [0.95, 0.83, 0.76], ...
                [1 1 1],...
                [0.8 0.9 0.96]};

CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));


%% Plot
%Establish grid: Should be 4x2

%Position of the plots organized in the same order as Files1
plotInds = [1, 3;... 
            5, 7;...
            2, 4;...
            6, 8];

xlims = [[10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3]];

xlabels = ["[MeAsp] \muM", "[MeAsp] \muM", "[Ser] \muM", "[Ser] \muM"];

figure()
for ii = 1:length(Files1)
    FilesTemp = Files1{ii};
    colr = colrs{ii};
    lnstyle = lnstyles{ii};
    lnwdth = lnwdths{ii};
    showMarker = showMarkers{ii};
    xlimits = xlims(ii, :);

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
        plot(LplotPDF, y ./ max(y), 'color', colr(jj), 'Linewidth', lnwdth(jj),...
            'linestyle', lnstyle(jj))
%         ylabel(['normalized', newline, 'pdf'])
        xlim([log(xlimits(1)), log(xlimits(2))])
        set(gca, 'XTick', [])
        if plotInds(ii, 2) /2 == round(plotInds(ii, 2)/2)
            set(gca, 'YTick', [])
        else
            ylabel("pdf")
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
        if showMarker(jj)
            shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(jj))
            errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr,...
                CDFPlotData.CDFPointPosErr, 'o', 'color', colr(jj));
        end
        C(jj) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), ...
            'color', colr(jj), 'Linewidth', lnwdth(jj), 'linestyle', lnstyle(jj));
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

set(gcf, 'Position', [1011.4,137.8,422.4000000000001,478.4])

saveas(gcf, './Outputs/Figure2.svg')
