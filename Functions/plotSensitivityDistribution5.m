function plotSensitivityDistribution(Files1, legend1, legend2, Files2, colr1, colr2, plotInd)
%% Declare fit functions
CDFFunc = @(p, L) logncdf(L, p(1), p(2));
% PDFFunc = @(p, L) lognpdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));
CV = []; %To store measurements of CV for each experiment
M = []; %To store measurements of mean K1/2 for each experiment
    
%% Load and plot CDF for all Data
figure()

%% MeAsp Dose response curves with different inhibitors
colr = colr1;
xlims = [10^-4.5, 10^3];


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

    %Plot CDF In Bottom left
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
        lgnd = legend(C, legend1, 'Position', [0.136518203933183,0.922560607161264,0.341541750591879,0.085227270525965]);
%         lgnd = legend(C, legend1, 'Location', 'Northwest'); 
        lgnd.FontSize = 16;
        lgnd.Orientation = 'horizontal'
        title(lgnd, "Foreground")
        set(lgnd,'color','none');
        legend boxoff
    end
    xlim(xlims)
    xticks([10^-4, 10^-2  10^0, 10^2, 10^4])
%     xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
    
    %Plot PDF in bottom left
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
%     shadeLineError(LplotPDF, curvePosErr2./max(y), curveNegErr2./max(y), colr(i))
%     plot(LplotPDF, y./max(y), 'color', colr(i), 'Linewidth', 1.8)

    shadeLineError(LplotPDF, curvePosErr2, curveNegErr2, colr(i))
    plot(LplotPDF, y, 'color', colr(i), 'Linewidth', 1.8)

    if i == length(Files1)
%         ylabel('p(K_{1/2})')
        ylabel('pdf')
        xlabel('log K_{1/2}')
%         sp.Position = sp.Position + [0 -0.07 0 -0.2];
    end
%     set(gca,'XTickLabel',[]);
    xlim([log(xlims)])
    
    %Calculate mean, variance, and CV from sampled lognormal parameters
    m = exp(rnd(:, 1) + rnd(:, 2).^2 ./ 2); %mean K1/2
    v = exp(2*rnd(:, 1) + rnd(:,2).^2).*(exp(rnd(:, 2).^2) - 1);
    CV = [CV, sqrt(v)./m];
    M = [M, m];
end

%% MeAsp Dose response curves with different inhibitors
colr = colr2;
xlims = [10^-4.5, 10^3];
for i = plotInd
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

    %Plot CDF In Bottom left
    subplot(2,2,2)
    hold on
    %Function to plot error of the curve as a fill using MCMC samples
    shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(i))

    errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr, CDFPlotData.CDFPointPosErr, 'o', 'color', colr(i));
    C(i) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), 'color', colr(i), 'Linewidth', 1.8);
    if i == length(Files2)
        xlabel('[MeAsp] \muM')
        ylabel('P(K_{1/2} < [L])')
        set(gca, 'xscale', 'log')
        
        %Handle legend position and font
        lgnd = legend(C, legend1, 'Position', [0.578392978010067,0.92298015068706,0.345238089883998,0.084388183474373]);
%         lgnd = legend(C, legend1, 'Location', 'Northwest'); 
        lgnd.FontSize = 16;
        lgnd.Orientation = 'horizontal'
        title(lgnd, "Background")
        set(lgnd,'color','none');
        legend boxoff
    end
    xlim(xlims)
    xticks([10^-4, 10^-2  10^0, 10^2, 10^4])
        
    %Plot PDF
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
%     shadeLineError(LplotPDF, curvePosErr2./max(y), curveNegErr2./max(y), colr(i))
%     plot(LplotPDF, y./max(y), 'color', colr(i), 'Linewidth', 1.8)

    shadeLineError(LplotPDF, curvePosErr2, curveNegErr2, colr(i))
    plot(LplotPDF, y, 'color', colr(i), 'Linewidth', 1.8)

%     if i == length(Files2)
%         ylabel('p(K_{1/2})')
        ylabel('pdf')
        xlabel('log K_{1/2}')
%         sp.Position = sp.Position + [0 -0.07 0 -0.2];
%     end
%     set(gca,'XTickLabel',[]);
    xlim(log(xlims))

    
    %Calculate mean, variance, and CV from sampled lognormal parameters
    m = exp(rnd(:, 1) + rnd(:, 2).^2 ./ 2); %mean K1/2
    v = exp(2*rnd(:, 1) + rnd(:,2).^2).*(exp(rnd(:, 2).^2) - 1);
    CV = [CV, sqrt(v)./m];
    M = [M, m];
end
JMAxes
% set(gcf, 'Position', [50, 50, 800, 600])
% set(gcf, 'Position', [74,43,1152,916])
set(gcf, 'Position', [86,38,924,711])

savefig(gcf, ['./Outputs/Figure3.fig'])
saveas(gcf, ['./Outputs/Figure3.png'])
saveas(gcf, ['./Outputs/Figure3.svg'])

%% Theoretical curves in the bottom left position

% %Generate distribution of Ntar and Ntsr
% muN = 9;
% sigmaN = 1;
% ncells = 100000;
% N1 = (lognrnd(2.0177, 0.3869, ncells, 1));
% N2 = (lognrnd(3.2322, 0.4174, ncells, 1));
% 
% %Generate distribution of K1/2
% Ki1 = 18; %uM
% Ki2 = 60; %uM
% alpha = 2; %kT
% m0 = 0.5;
% mstar = 0.4450;
% 
% 
% %Distribution of K1/2 for Asp
% KhalfAsp = @(Asp) Ki1*((exp((N1)*alpha*(m0-mstar)) +2).^(1./N1)* ...
%     (1+Asp./Ki1) -1)
% 
% KhalfAsp2 = @(Asp) Ki2*((exp((N1)*alpha*(m0-mstar)) +2).^(1./N1)* ...
%     (1+Asp./Ki2) -1)
% JMAxes
% subplot(2,2, 3)
% hold on
% %histogram(log(KhalfAsp(0)), 'Normalization',"probability")
% [f1, xi1] = ksdensity(log(KhalfAsp(0)));
% plot(xi1, f1, 'LineWidth', 6, 'color', colr(1))
% plot(xi1, f1, 'LineWidth', 2, 'color', colr(2))
% 
% %histogram(log(KhalfAsp2(0)),'Normalization',"probability")
% [f2, xi2] = ksdensity(log(KhalfAsp2(0)));
% plot(xi2, f2, 'Linewidth', 3, 'color', colr(3))
% xlabel("log(K_{1/2}^{Asp})")
% lgnd = legend(["0 Background", "Uncompetetive", "Competetive"]);
% % title(lgnd, "Background")
%         set(lgnd,'color','none');
%         legend boxoff
%         lgnd.Position = [0.142916664217079,0.452555552980635,0.204750003361702,0.116333335908254]
% set(gca, 'FontSize', 14)
% set(gcf, 'Position', [50, 50, 800, 600])
% set(gca,'YTickLabel',[], 'XTickLabel', []);
% xlabel("Log[L]")
% ylabel('P(K_{1/2} = [L])')


%% Plot CV measurements and confidence intervals for each condition
figure()
avgCV = mean(CV, 1);
CVposErr = prctile(CV, 97.5, 1);
CVnegErr = prctile(CV, 2.5, 1);

subplot(2, 2, 1)
errorbar(1:length(avgCV),avgCV, avgCV - CVnegErr, CVposErr - avgCV, 'o')
xlim([0.5, 3.5])
xticks([1, 2, 3, 4, 5, 6])
xticklabels([legend1, legend2])
hold on
ylim([0,3.5])
ylabel("CV(K_{1/2})")

subplot(2, 2, 2)
errorbar(1:length(avgCV),avgCV, avgCV - CVnegErr, CVposErr - avgCV, 'o')
xlim([3.5, 6.5])
xticks([1, 2, 3, 4, 5, 6])
xticklabels([legend1, legend2])
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
xticklabels([legend1, legend2])

subplot(2, 2, 4)
avgM = mean(M, 1);
MposErr = prctile(M, 97.5, 1);
MnegErr = prctile(M, 2.5, 1); 
errorbar(1:length(avgM),avgM, avgM - MnegErr, MposErr - avgM, 'o')
set(gca, 'yscale', 'log')
xlim([3.5, 6.5])
xticks([1, 2, 3, 4, 5, 6])
xticklabels([legend1, legend2])

JMAxes
set(gcf, 'Position', [489,289,572,474])
end