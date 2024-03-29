function plotSensitivityDistribution7(Files1, legend1, colr, plotInd)
%% Declare fit functions
CDFFunc = @(p, L) logncdf(L, p(1), p(2));
% PDFFunc = @(p, L) lognpdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));
CV = []; %To store measurements of CV for each experiment
M = []; %To store measurements of mean K1/2 for each experiment
    
%% Load and plot CDF for all Data
figure()
plotNum = 1;

for i = plotInd
    %Load Dataset
    load([convertStringsToChars(Files1(i)), 'plotData.mat'])
            Lplot = CDFPlotData.Lplot;
        %Calculate error bar on CDF curve
        %Generate individual curves for each of the MCMC Samples and calculate
        %95-CI
        rnd = CDFPlotData.MCMCSamples;
        curveSamples = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples(j, :) = CDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr = prctile(curveSamples, 97.5);
        curveNegErr = prctile(curveSamples, 2.5);

    %Plot CDF In Top left
    subplot(2,2,1)
    hold on
            %Function to plot error of the curve as a fill using MCMC samples
    shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(i))

    errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr, CDFPlotData.CDFPointPosErr, 'o', 'color', colr(i));
    C(plotNum) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), 'color', colr(i), 'Linewidth', 1.8);
    if i == plotInd(end)
        xlabel('[meAsp] \muM')
        ylabel('P(K_{1/2} < [L])')
        set(gca, 'xscale', 'log')
        
        %Handle legend position and font
%         lgnd = legend(C, legend1, 'Position', [.55, .80, .2, .1]);
        lgnd = legend(C, legend1, 'Position', [0.136518203933183,0.922560607161264,0.341541750591879,0.085227270525965]); 
        lgnd.FontSize = 16;
        lgnd.Orientation = 'horizontal';
        title(lgnd, "Background")
        set(lgnd,'color','none');
        legend boxoff
    end
    
        curveSamples2 = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples2(j, :) = PDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr2 = prctile(curveSamples2, 97.5);
        curveNegErr2 = prctile(curveSamples2, 2.5);
    xlim([10^0, 10^3])



    %Plot PDF in top left
    LplotPDF = [log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot))];
    curveSamples2 = zeros(size(rnd,1), length(LplotPDF));
    for j = 1:size(rnd, 1)
        curveSamples2(j, :) = PDFFunc(rnd(j, :), LplotPDF);
    end
    curvePosErr2 = prctile(curveSamples2, 97.5);
    curveNegErr2 = prctile(curveSamples2, 2.5);
    
    subplot(2, 2, 3)
    hold on
    y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
    shadeLineError(LplotPDF, curvePosErr2, curveNegErr2, colr(i))
    plot(LplotPDF, y, 'color', colr(i), 'Linewidth', 1.8)
    xlabel('log K_{1/2}')
    ylabel('pdf')
    xlim([log(10^0), log(10^3)])

    
    
    %Calculate mean, variance, and CV from sampled lognormal parameters
    m = exp(rnd(:, 1) + rnd(:, 2).^2 ./ 2); %mean K1/2
    v = exp(2*rnd(:, 1) + rnd(:,2).^2).*(exp(rnd(:, 2).^2) - 1);
    CV = [CV, sqrt(v)./m];
    M = [M, m];
    
    plotNum = plotNum + 1;
end
JMAxes
% set(gcf, 'Position', [50, 50, 800, 600])
% set(gcf, 'Position', [74,43,1152,916])
set(gcf, 'Position', [86,38,924,711])



end