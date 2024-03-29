function plotFRETTimeSeries(EfretData, dataparms, parms, OutputDest)
    %% Reorganize data to make it amenable to plotting time-series
    plotData = getPlotFriendlyData(EfretData);
    avgASeries = mean(plotData.aSeries, 1);
    avgTSeries = mean(plotData.tSeries, 1);
    switchTimes = find(plotData.sChanges);
    
    %% Plot population-average time-series
    avgTSeries = avgTSeries - min(avgTSeries);
    for i = 1:length(avgASeries)
       if i / 20 == round(i/20)
           avgTSeries(i) = nan;
       end
    end

    figure(); hold on
    % for i = 1:size(plotData.aSeries, 1)
    %     plot(avgTSeries, smooth(plotData.aSeries(i, :)), 'Color', [.8 .8 .8])
    % end
    subplot(2, 1, 1)
    plot(avgTSeries, avgASeries, 'Color', parms.lineColor, 'Linewidth', 3)
    xlabel('Time (s)')
    ylabel('\langle a \rangle')
    set(gca, {'XColor', 'YColor'}, {[0.4 0.4 0.4], [0.4 0.4 0.4]});
    set(gca, 'FontSize', 20)
    
    %Color stimulus regions
    for i = 1:2:length(switchTimes)
        if i < length(switchTimes)
            x = avgTSeries(transpose(switchTimes(i:i+1) - [0;2])); %Had to subtract 2 because of nans
        else
            x = [avgTSeries(switchTimes(i)), max(avgTSeries)];
        end
        
        y = get(gca,'YLim');
        nx = [x(1) x(2) x(2) x(1)];
        ny = [y(1) y(1) y(2) y(2)];
        hold on
        patch(nx,ny, 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end

    %Plot in terms of frames
    frames = 1:length(avgASeries);
    subplot(2, 1, 2)
    plot(frames, avgASeries, 'Color', parms.lineColor, 'Linewidth', 3)
    xlabel('Frames')
    ylabel('\langle a \rangle')
    set(gca, {'XColor', 'YColor'}, {[0.4 0.4 0.4], [0.4 0.4 0.4]});
    set(gca, 'FontSize', 20)
    
  %Color stimulus regions
    for i = 1:2:length(switchTimes)
        if i < length(switchTimes)
            x = frames(transpose(switchTimes(i:i+1) - [0;1])); %subtract 1 because no nans
        else
            x = [frames(switchTimes(i)), max(frames)]; 
        end
        
        y = get(gca,'YLim');
        nx = [x(1) x(2) x(2) x(1)];
        ny = [y(1) y(1) y(2) y(2)];
        hold on
        patch(nx,ny, 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    end

    set(gcf, 'Position', [200,100, 1000, 600])
    
    %Save figure
    savefig(gcf, [OutputDest, 'PopAvgTimeSeries.fig'])
    saveas(gcf, [OutputDest, 'PopAvgTimeSeries.png'])

    
    %% Plot single-cell traces
    mkdir([OutputDest, 'SingleCellTraces/'])
    
    nPlots = 10;
    for i = 1:size(plotData.aSeries, 1)
        if mod(i, nPlots) == 1
            figure('Visible', 'on') 
        end
        %Plot single-cell trace
        subplot(nPlots./2, 2, mod(i-1, nPlots)+1)
        hold on
        plot(plotData.aSeries(i, :), 'o', 'color', [0.4, 0.4, 0.4]);
        plot(smooth(plotData.aSeries(i, :)), 'r', 'Linewidth', 2);
        ylabel("a")
        ylim([-0.2, .8])
         %Color stimulus regions
        for k = 1:2:length(switchTimes)
            if k < length(switchTimes)
                x = frames(transpose(switchTimes(k:k+1) - [0;1])); %subtract 1 because no nans
            else
                x = [frames(switchTimes(k)), max(frames)]; 
            end

            y = get(gca,'YLim');
            nx = [x(1) x(2) x(2) x(1)];
            ny = [y(1) y(1) y(2) y(2)];
            hold on
            p(k) = patch(nx,ny, 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            uistack(p(k), 'bottom')
        end
        
        if ~(mod(i, nPlots) == 0 | i == size(plotData.aSeries, 1) | mod(i, nPlots) == nPlots-1)
           set(gca, 'xtick', []) 
        end
        
        if mod(i, nPlots)==0 | mod(i, nPlots) == nPlots-1 | i == size(plotData.aSeries, 1)
            xlabel('Frames')
        end
        
        if mod(i, nPlots) == 0 | i == size(plotData.aSeries, 1)
            set(gcf, 'Position', [100,25, 1200, 750])
            savefig(gcf, [OutputDest, 'SingleCellTraces/upto', num2str(i), '.fig'])
            saveas(gcf, [OutputDest, 'SingleCellTraces/upto', num2str(i), '.png'])
        end
    end
  


end