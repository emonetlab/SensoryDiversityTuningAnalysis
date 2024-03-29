function plotSensitivityDistribution(Files1, legend1, Files2, colr, plotInd)
%% Declare fit functions
CDFFunc = @(p, L) logncdf(L, p(1), p(2));
% PDFFunc = @(p, L) lognpdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));
CV = []; %To store measurements of CV for each experiment
M = []; %To store measurements of mean K1/2 for each experiment
    
%% Load and plot CDF for all Data
figure()

%% MeAsp K1/2 distributions
for i = plotInd
    %Load Dataset
    load([convertStringsToChars(Files1(i)), 'plotData.mat'])
        Lplot = CDFPlotData.Lplot;
        %Calculate error bar on CDF curve
        %Generate individual curves for each of the MCMC Samples and calculate
        %95-CI
        rnd = CDFPlotData.MCMCSamples; %Dimensions: [mu, sigma, A (not used)]
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
    xlabel('[MeAsp] \muM')
    ylabel('P(K_{1/2} < [L])')
    set(gca, 'xscale', 'log')
    xlim([10^-3, 10^3])
    
    %Plot PDF in top left
    LplotPDF = [log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot))];
    curveSamples2 = zeros(size(rnd,1), length(LplotPDF));
    for j = 1:size(rnd, 1)
        curveSamples2(j, :) = PDFFunc(rnd(j, :), LplotPDF);
    end
    curvePosErr2 = prctile(curveSamples2, 97.5);
    curveNegErr2 = prctile(curveSamples2, 2.5);
    
    sp = subplot(2, 2, 3);
    hold on
    y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
%     shadeLineError(CDFPlotData.Lplot, curvePosErr2./max(y), curveNegErr2./max(y), colr(i))
%     plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)

    shadeLineError(LplotPDF, curvePosErr2, curveNegErr2, colr(i))
    plot(LplotPDF, y, 'color', colr(i), 'Linewidth', 1.8)
    ylabel('pdf')
    xlim([log(10^-3), log(10^3)])
    xlabel('log K_{1/2}')
    
    
    %Calculate mean, variance, and CV from sampled lognormal parameters
    m = exp(rnd(:, 1) + rnd(:, 2).^2 ./ 2); %mean K1/2
    v = exp(2*rnd(:, 1) + rnd(:,2).^2).*(exp(rnd(:, 2).^2) - 1);
    CV = [CV, sqrt(v)./m];
    M = [M, m];
end

%% Ser K1/2 distributions
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

    %Plot CDF In top right
    subplot(2,2,2)
    hold on
            %Function to plot error of the curve as a fill using MCMC samples
    shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(i))

    errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr, CDFPlotData.CDFPointPosErr, 'o', 'color', colr(i));
    C(i) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), 'color', colr(i), 'Linewidth', 1.8);
    if i == length(Files2)
        xlabel('[Ser] \muM')
        ylabel('P(K_{1/2} < [L])')
        set(gca, 'xscale', 'log')
        
        %Handle legend position and font
%         lgnd = legend(C, legend1, 'Position', [.55, .80, .2, .1]);
        lgnd = legend(C, legend1, 'Location', 'Northwest'); 
        lgnd.FontSize = 16;
        lgnd.Position = [0.150095535541953,0.923657383727201,0.341991336734006,0.084388183474373];
        lgnd.Orientation = 'horizontal';
        title(lgnd, "Background")
        set(lgnd,'color','none');
        legend boxoff
    end
    xlim([10^-4, 10^2])
    
    %Plot PDF in bottom right
    LplotPDF = [log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot))];
    curveSamples2 = zeros(size(rnd,1), length(LplotPDF));
    for j = 1:size(rnd, 1)
        curveSamples2(j, :) = PDFFunc(rnd(j, :), LplotPDF);
    end
    curvePosErr2 = prctile(curveSamples2, 97.5);
    curveNegErr2 = prctile(curveSamples2, 2.5);

    sp = subplot(2, 2, 4);
    hold on
    y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
%     shadeLineError(CDFPlotData.Lplot, curvePosErr2./max(y), curveNegErr2./max(y), colr(i))
%     plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)

    shadeLineError(LplotPDF, curvePosErr2, curveNegErr2./max(y), colr(i))
    plot(LplotPDF, y, 'color', colr(i), 'Linewidth', 1.8)

    if i == length(Files2)
%         xlabel('[Ser] \muM')
%         ylabel('p(log K_{1/2})')
          ylabel('pdf')
%         set(gca, 'xscale', 'log')
%         sp.Position = sp.Position + [0 -0.07 0 -0.2];
    end
    xlim([log(10^-4), log(10^2)])
%     set(gca,'XTickLabel',[]);
    xlabel('log K_{1/2}')
    
    %Calculate mean, variance, and CV from sampled lognormal parameters
    m = exp(rnd(:, 1) + rnd(:, 2).^2 ./ 2); %mean K1/2
    v = exp(2*rnd(:, 1) + rnd(:,2).^2).*(exp(rnd(:, 2).^2) - 1);
    CV = [CV, sqrt(v)./m];
    M = [M, m];
end

% set(gcf, 'Position', [50, 50, 800, 600])
% set(gcf, 'Position', [74,43,1152,916])
set(gcf, 'Position', [86,38,924,711])
JMAxes
savefig(gcf, ['./Outputs/Figure2.fig'])
saveas(gcf, ['./Outputs/Figure2.png'])
saveas(gcf, ['./Outputs/Figure2.svg'])


%% Plot insets of 0 background cases.
figure()
%MeAsp
for i = [3, 1]
    %Load Dataset
    load([convertStringsToChars(Files1(i)), 'plotData.mat'])
    Lplot = CDFPlotData.Lplot;
    %Calculate error bar on CDF curve
    %Generate individual curves for each of the MCMC Samples and calculate
    %95-CI
    rnd = CDFPlotData.MCMCSamples; %Dimensions: [mu, sigma, A (not used)]
    curveSamples = zeros(size(rnd,1), length(Lplot));
    for j = 1:size(rnd, 1)
        curveSamples(j, :) = CDFFunc(rnd(j, :), Lplot);
    end
    curvePosErr = prctile(curveSamples, 97.5);
    curveNegErr = prctile(curveSamples, 2.5);
    
    %Plot PDF in top left
    LplotPDF = [log(10^-3):0.01:log(10^3)];
    curveSamples2 = zeros(size(rnd,1), length(LplotPDF));
    for j = 1:size(rnd, 1)
        curveSamples2(j, :) = PDFFunc(rnd(j, :), LplotPDF);
    end
    curvePosErr2 = prctile(curveSamples2, 97.5);
    curveNegErr2 = prctile(curveSamples2, 2.5);
    
    subplot(1, 2, 1);
    hold on
    y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
    shadeLineError(LplotPDF, curvePosErr2, curveNegErr2, colr(i))
    plot(LplotPDF, y, 'color', colr(i), 'Linewidth', 1.8)
%     xlabel("log K_{1/2}")
    xlim([min(LplotPDF), max(LplotPDF)])
end

%Serine
for i = [1, 2]
    %Load Dataset
    load([convertStringsToChars(Files2(i)), 'plotData.mat'])
    Lplot = CDFPlotData.Lplot;
    %Calculate error bar on CDF curve
    %Generate individual curves for each of the MCMC Samples and calculate
    %95-CI
    rnd = CDFPlotData.MCMCSamples; %Dimensions: [mu, sigma, A (not used)]
    curveSamples = zeros(size(rnd,1), length(Lplot));
    for j = 1:size(rnd, 1)
        curveSamples(j, :) = CDFFunc(rnd(j, :), Lplot);
    end
    curvePosErr = prctile(curveSamples, 97.5);
    curveNegErr = prctile(curveSamples, 2.5);
    
    %Plot PDF in top left
    LplotPDF = [log(10^-3):0.01:log(10^3)];
    curveSamples2 = zeros(size(rnd,1), length(LplotPDF));
    for j = 1:size(rnd, 1)
        curveSamples2(j, :) = PDFFunc(rnd(j, :), LplotPDF);
    end
    curvePosErr2 = prctile(curveSamples2, 97.5);
    curveNegErr2 = prctile(curveSamples2, 2.5);
    
    subplot(1, 2, 2);
    hold on
    y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
    shadeLineError(LplotPDF, curvePosErr2, curveNegErr2, colr(i))
    plot(LplotPDF, y, 'color', colr(i), 'Linewidth', 1.8)
%     xlabel("log K_{1/2}")
    xlim([min(LplotPDF), 1])
end
JMAxes
set(gcf, 'Position', [600,343,560,126.8])
% set(gcf, 'Position', [69,50,1128,902])


%% Plot CV measurements and confidence intervals for each condition
figure()
avgCV = mean(CV, 1);
CVposErr = prctile(CV, 97.5, 1);
CVnegErr = prctile(CV, 2.5, 1);

subplot(2, 2, 1)
errorbar(1:length(avgCV),avgCV, avgCV - CVnegErr, CVposErr - avgCV, 'o')
xlim([0.5, 3.5])
xticks([1, 2, 3, 4, 5, 6])
xticklabels([legend1, legend1])
hold on
ylim([0,3.5])
ylabel("CV(K_{1/2})")

subplot(2, 2, 2)
errorbar(1:length(avgCV),avgCV, avgCV - CVnegErr, CVposErr - avgCV, 'o')
xlim([3.5, 6.5])
xticks([1, 2, 3, 4, 5, 6])
xticklabels([legend1, legend1])
ylim([0,3.5])


subplot(2, 2, 3)
avgM = mean(M, 1);
MposErr = prctile(M, 97.5, 1);
MnegErr = prctile(M, 2.5, 1); 
errorbar(1:length(avgM),avgM, avgM - MnegErr, MposErr - avgM, 'o')
set(gca, 'yscale', 'log')
xlim([.5, 3.5])
ylabel("\langleK_{1/2}\rangle")
xticks([1, 2, 3, 4, 5, 6])
xticklabels([legend1, legend1])

subplot(2, 2, 4)
avgM = mean(M, 1);
MposErr = prctile(M, 97.5, 1);
MnegErr = prctile(M, 2.5, 1); 
errorbar(1:length(avgM),avgM, avgM - MnegErr, MposErr - avgM, 'o')
set(gca, 'yscale', 'log')
xlim([3.5, 6.5])
xticks([1, 2, 3, 4, 5, 6])
xticklabels([legend1, legend1])

JMAxes
set(gcf, 'Position', [489,289,572,474])



end