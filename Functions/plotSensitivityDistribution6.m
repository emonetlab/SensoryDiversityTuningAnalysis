function plotSensitivityDistribution(Files1, legend1, legend2, Files2, colr)
%% Declare fit functions
CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) lognpdf(L, p(1), p(2));
CV = []; %To store measurements of CV for each experiment
M = []; %To store measurements of mean K1/2 for each experiment
    
%% Load and plot CDF for all Data
figure()

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

    errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr, CDFPlotData.CDFPointPosErr, 'o', 'color', colr(i));
    C(i) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), 'color', colr(i), 'Linewidth', 1.8);
    if i == length(Files1)
        xlabel('[L] \muM')
        ylabel('P(K_{1/2} < [L])')
        set(gca, 'xscale', 'log')
        
        %Handle legend position and font
%         lgnd = legend(C, legend1, 'Position', [.55, .80, .2, .1]);
        lgnd = legend(C, legend1, 'Location', 'Northwest'); 
        lgnd.FontSize = 16;
%         title(lgnd, "Background")
        set(lgnd,'color','none');
        legend boxoff
    end
    
        curveSamples2 = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples2(j, :) = PDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr2 = prctile(curveSamples2, 97.5);
        curveNegErr2 = prctile(curveSamples2, 2.5);
    xlim([10^-4, 10^4])

    subplot(2, 2, 3)
    hold on
    y = PDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot);
    shadeLineError(CDFPlotData.Lplot, curvePosErr2./max(y), curveNegErr2./max(y), colr(i))
    plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)

    plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)
    if i == length(Files1)
        xlabel('[L] \muM')
        ylabel('P(K_{1/2} = [L])')
        set(gca, 'xscale', 'log')
    end
    xlim([10^-4, 10^4])

    
    %Calculate mean, variance, and CV from sampled lognormal parameters
    m = exp(rnd(:, 1) + rnd(:, 2).^2 ./ 2); %mean K1/2
    v = exp(2*rnd(:, 1) + rnd(:,2).^2).*(exp(rnd(:, 2).^2) - 1);
    CV = [CV, sqrt(v)./m];
    M = [M, m];
end


%% Dose response curves for different Tar-binding ligands
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
    C(i) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), 'color', colr(i), 'Linewidth', 1.8);
    if i == length(Files2)
        xlabel('[L] \muM')
        ylabel('P(K_{1/2} < [L])')
        set(gca, 'xscale', 'log')
        
        %Handle legend position and font
        lgnd = legend(C, legend2, 'Position', [0.84,0.77,0.14,0.13]);
%         lgnd = legend(C, legend2, 'Location', 'Northwest'); 
        lgnd.FontSize = 16;
%         title(lgnd, "Background")
        set(lgnd,'color','none');
        legend boxoff
    end
    
        curveSamples2 = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples2(j, :) = PDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr2 = prctile(curveSamples2, 97.5);
        curveNegErr2 = prctile(curveSamples2, 2.5);
    xlim([10^-4, 10^4])

    subplot(2, 2, 4)
    hold on
    y = PDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot);
    shadeLineError(CDFPlotData.Lplot, curvePosErr2./max(y), curveNegErr2./max(y), colr(i))
    plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)

    plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)
    if i == length(Files2)
        xlabel('[L] \muM')
        ylabel('P(K_{1/2} = [L])')
        set(gca, 'xscale', 'log')
    end
    xlim([10^-4, 10^4])

    
    %Calculate mean, variance, and CV from sampled lognormal parameters
    m = exp(rnd(:, 1) + rnd(:, 2).^2 ./ 2); %mean K1/2
    v = exp(2*rnd(:, 1) + rnd(:,2).^2).*(exp(rnd(:, 2).^2) - 1);
    CV = [CV, sqrt(v)./m];
    M = [M, m];
end
JMAxes
set(gcf, 'Position', [50, 50, 800, 600])



%% Plot CV measurements and confidence intervals for each condition
figure()
avgCV = mean(CV, 1);
CVposErr = prctile(CV, 97.5, 1);
CVnegErr = prctile(CV, 2.5, 1);
subplot(2, 1, 1)
errorbar(1:length(avgCV),avgCV, avgCV - CVnegErr, CVposErr - avgCV, 'o')
xlim([0, length(legend1)+0.5])
% set(gca, 'yscale', 'log')
title("CV")

subplot(2, 1, 2)
avgM = mean(M, 1);
MposErr = prctile(M, 97.5, 1);
MnegErr = prctile(M, 2.5, 1); 
errorbar(1:length(avgM),avgM, avgM - MnegErr, MposErr - avgM, 'o')
set(gca, 'yscale', 'log')
xlim([0, length(legend1)+0.5])
title("Mean K_{1/2}")
JMAxes



end