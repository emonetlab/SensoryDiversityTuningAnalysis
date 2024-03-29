function plotSensitivityDistribution2(Files1, legend1)
%% Declare fit functions
CDFFunc = @(p, L) p(3)*logncdf(L, p(1), p(2));
PDFFunc = @(p, L) lognpdf(L, p(1), p(2));

%% Color Palette To Color Each Background
colr = ["#058ED9", "#FF934F", "#CC2D35", "#2D3142", "#848FA2"];
    
%% Load and plot CDF for all Data
figure()

%MeAsp dose response curves
for i = 1:length(Files1)
    %Load Dataset
    load([convertStringsToChars(Files1(i)), 'plotData.mat'])
    
    %Plot CDF In Top left
    subplot(2,1,1)
    hold on
    errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointStd, 'o', 'color', colr(i));
    C(i) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), 'color', colr(i), 'Linewidth', 1.8);
    if i == length(Files1)
        xlabel('[MeAsp] \muM')
        ylabel('P(K_{1/2} < [L])')
        set(gca, 'xscale', 'log')
    end
    lgnd = legend(C, legend1, 'Location', 'Northwest'); 
        lgnd.FontSize = 16;
        set(lgnd,'color','none');
        legend boxoff
    
    subplot(2, 1, 2)
    hold on
    y = PDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot);
    plot(CDFPlotData.Lplot, y./max(y), 'color', colr(i), 'Linewidth', 1.8)
    if i == length(Files1)
        xlabel('[MeAsp] \muM')
        ylabel('P(K_{1/2} = [L])')
        set(gca, 'xscale', 'log')
    end
    
end


set(gcf, 'Position', [50, 50, 500, 600])


end