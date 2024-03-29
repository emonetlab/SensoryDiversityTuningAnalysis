function plotFretRandDistribution(Files1, legend1, Files2)
    %% Color Palette To Color Each Background
    colr = ["#058ED9", "#FF934F", "#CC2D35"];

    %% Load and plot CDF for all Data
    figure()
    
    for i = 1:length(Files1)
        %Load Dataset
        load([convertStringsToChars(Files1(i)), 'plotData.mat'])
        subplot(2, 2, 1)
        title("MeAsp Dose Response Curves")
        hold on
        C(i) = histogram(FretRangeData.diff_all, 'FaceColor', colr(i),'Normalization','probability');
        if i == length(Files1)
            lgnd = legend(C, legend1, 'Location', 'best'); 
            set(lgnd,'color','none');
            legend boxoff
            xlabel("FretMax - FretMin")

        end

        
        subplot(2, 2, 3)
        hold on
        histogram(FretRangeData.diffNorm_all, 'FaceColor', colr(i), 'Normalization','probability')
        xlabel("(FretMax - FretMin)/FretMin")

    end
    
        for i = 1:length(Files2)
        %Load Dataset
        load([convertStringsToChars(Files2(i)), 'plotData.mat'])
        subplot(2, 2, 2)
        title("Ser Dose Response Curves")
        hold on
        C(i) = histogram(FretRangeData.diff_all, 'FaceColor', colr(i),'Normalization','probability');
        if i == length(Files1)
            xlabel("FretMax - FretMin")

        end

        
        subplot(2, 2, 4)
        hold on
        histogram(FretRangeData.diffNorm_all, 'FaceColor', colr(i), 'Normalization','probability')
        xlabel("(FretMax - FretMin)/FretMin")
    end
    set(gcf, 'Position', [50, 50, 800, 600])


end