function plotData = getPlotFriendlyData(EfretData)   
    %Get long time-series for each cell
    for i = 1:length(EfretData)
        a = transpose(EfretData(i).a);
        t = transpose(EfretData(i).t);
        plotData.aSeries(i, :) = a(:);
        plotData.tSeries(i, :) = t(:);
        
        a_sat1 = transpose(EfretData(i).a_min_max(1, :));
        t_sat1 = transpose(EfretData(i).t_min_max(1,:));
        a_sat2 = transpose(EfretData(i).a_min_max(2, :));
        t_sat2 = transpose(EfretData(i).t_min_max(2,:));
        plotData.a_sat1Series(i, :) = a_sat1(:);
        plotData.t_sat1Series(i, :) = t_sat1(:);
        plotData.a_sat2Series(i, :) = a_sat2(:);
        plotData.t_sat2Series(i, :) = t_sat2(:);
    end
    
    %Get indices where stimulus changes, now updated to include minmax
    s = transpose(EfretData(1).s);
    plotData.sSeries = s(:);
    
    s_sat1 = transpose(EfretData(1).s_min_max(1, :));
    s_sat2 = transpose(EfretData(1).s_min_max(2, :));
    plotData.s_sat1Series = s_sat1(:);
    plotData.s_sat2Series = s_sat2(:);
    
    plotData.sChanges = zeros(length(plotData.sSeries), 1);
    for ii = 2:length(plotData.sSeries)
        if plotData.sSeries(ii) ~= plotData.sSeries(ii-1)
            plotData.sChanges(ii) = 1;
        end
    end
    
    plotData.s_sat1Changes = zeros(length(plotData.s_sat1Series), 1);
    for ii = 2:length(plotData.s_sat1Series)
        if plotData.s_sat1Series(ii) ~= plotData.s_sat1Series(ii-1)
            plotData.s_sat1Changes(ii) = 1;
        end
    end
    
    plotData.s_sat2Changes = zeros(length(plotData.s_sat2Series), 1);
    for ii = 2:length(plotData.s_sat2Series)
        if plotData.s_sat2Series(ii) ~= plotData.s_sat2Series(ii-1)
            plotData.s_sat2Changes(ii) = 1;
        end
    end
end