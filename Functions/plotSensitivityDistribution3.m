function plotSensitivityDistribution(Files1, legend1, Files2, colr)
%% Declare fit functions
CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) lognpdf(L, p(1), p(2));


    
%% Load and plot CDF for all Data
figure()

%MeAsp dose response curves
for i = 1:length(Files1)
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

    C(i) = errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr, CDFPlotData.CDFPointPosErr, 'o', 'color', colr(i));
    plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), 'color', colr(i), 'Linewidth', 1.8);
    if i == length(Files1)
        xlabel('[Maltose] \muM')
        ylabel('P(K_{1/2} < [L])')
        set(gca, 'xscale', 'log')
    end
            curveSamples2 = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples2(j, :) = PDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr2 = prctile(curveSamples2, 97.5);
        curveNegErr2 = prctile(curveSamples2, 2.5);

    subplot(2, 2, 3)
    hold on
    y = PDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot);
        shadeLineError(CDFPlotData.Lplot, curvePosErr2./max(y), curveNegErr2./max(y), colr(i))
    plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)

    plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)
    if i == length(Files1)
        xlabel('[Maltose] \muM')
        ylabel('P(K_{1/2} = [L])')
        set(gca, 'xscale', 'log')
    end
end

%Ser Dose-response curves
for i = 1:length(Files2)
    %Load Dataset
    load([convertStringsToChars(Files2(i)), 'plotData.mat'])
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
    subplot(2,2,2)
    hold on
            %Function to plot error of the curve as a fill using MCMC samples
        shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(i))

     errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr, CDFPlotData.CDFPointPosErr, 'o', 'color', colr(i));
    C(i) =plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), 'color', colr(i), 'Linewidth', 1.8);
    if i == length(Files2)
        xlabel('[MeAsp] \muM')
        ylabel('P(K_{1/2} < [L])')
        set(gca, 'xscale', 'log')
        
        %Handle legend position and font
%         lgnd = legend(C, legend1, 'Position', [.55, .80, .2, .1]);
        lgnd = legend(C, legend1, 'Location', 'Northwest'); 
        lgnd.FontSize = 16;
        set(lgnd,'color','none');
        legend boxoff
    end
    
        curveSamples2 = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples2(j, :) = PDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr2 = prctile(curveSamples2, 97.5);
        curveNegErr2 = prctile(curveSamples2, 2.5);

    subplot(2, 2, 4)
    hold on
    y = PDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot);
    shadeLineError(CDFPlotData.Lplot, curvePosErr2./max(y), curveNegErr2./max(y), colr(i))
    plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)

    plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)
    if i == length(Files1)
        xlabel('[MeAsp] \muM')
        ylabel('P(K_{1/2} = [L])')
        set(gca, 'xscale', 'log')
    end
end

set(gcf, 'Position', [50, 50, 800, 600])


end