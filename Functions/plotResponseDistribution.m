function plotResponseDistribution(EfretData, dataparms, parms, OutputDest, m)
%% Boxplot of all response levels
    %Reformat data to make easier for plotting
    plotData = getPlotFriendlyData(EfretData);
    
    %First, find all the stimulus concentrations
    stimlevels = [];
    for i = 1:length(plotData.sSeries)
        if ~ismember(plotData.sSeries(i), stimlevels)
            stimlevels = [stimlevels, plotData.sSeries(i)];
        end
    end
    %Remove the background level
    stimlevels(find(stimlevels == min(stimlevels))) = [];
    
    %Get all FRET values for each stimulus level
    figure();
    for i = 1:length(stimlevels)
        %Get values for each stimulus level
        fretVals = plotData.aSeries(:, plotData.sSeries == stimlevels(i));
        %Sort into columns based on how many stimuli were applied
        for j = 1:m
            temp = fretVals(:, (j-1) * size(fretVals, 2)./m + 1:j * size(fretVals, 2)./m);
            temp = mean(temp, 2); %If you want within-cell average
            sortedFretVals(:, j) = temp(:);
        end
        
        subplot(length(stimlevels), 1, i)
        boxplot(sortedFretVals)
        hold on
        yline(median(sortedFretVals(:)))
        title(['[L] = ', num2str(stimlevels(i))])
        ylabel('a')
        
        if i ~= length(stimlevels)
            set(gca,'XTick',[])
        else
            xlabel("Response Number")
        end
        
        set(gcf, 'Position', [100, 0, 600, 800])
        
    end
    
    savefig(gcf, [OutputDest, 'ResponseDistribution.fig'])
    saveas(gcf, [OutputDest, 'ResponseDistribution.png'])
    saveas(gcf, [OutputDest, 'ResponseDistribution.svg'])

%% Histogram of strongest response level, and pre-stimulus levels
    %First, find all the stimulus concentrations
    stimlevels = [];
    for i = 1:length(plotData.sSeries)
        if ~ismember(plotData.sSeries(i), stimlevels)
            stimlevels = [stimlevels, plotData.sSeries(i)];
        end
    end

    %Background kinase activity
    fretVals_bckgr = plotData.aSeries(:, plotData.sSeries == min(stimlevels));
    fretVals_max = plotData.aSeries(:, plotData.sSeries == max(stimlevels));

    figure()
    hold on
    histogram(fretVals_bckgr, 'FaceAlpha', 0.3, 'normalization', 'pdf')
    histogram(fretVals_max, 'FaceAlpha', 0.3, 'normalization', 'pdf')
    legend(["Background", "Max Stimulus"])
    xlabel("a")
    ylabel("pdf")
    savefig(gcf, [OutputDest, 'ResponseHistogram.fig'])
    saveas(gcf, [OutputDest, 'ResponseHistogram.png'])


end