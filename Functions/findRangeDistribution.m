function FretRange = findRangeDistribution(EfretData, OutputDest)
    %% Calculate range of fret values for each cell
    nCells = length(EfretData);
    diff_all = zeros(nCells, 1);
    diffNorm_all = zeros(nCells, 1);
    
    for i = 1:nCells
        %Get min and max fret values
        s = EfretData(i).s_min_max;
        Efret = EfretData(i).FRET_min_max;
        smin = (s == max(s(:)));
        smax = (s == min(s(:)));
        EfretMin = Efret(smin);
        EfretMax = Efret(smax);
        
        %Calculate difference and normalized difference
        diff_all(i) = median(EfretMax) - median(EfretMin);
        diffNorm_all(i) = (median(EfretMax) - median(EfretMin))/median(EfretMin);
%         diff_all(i) = mean(EfretMax) - mean(EfretMin);
%         diffNorm_all(i) = (mean(EfretMax) - mean(EfretMin))/mean(EfretMin);
    end
    
    %% Plot distributions of these ranges
    figure()
    subplot(2, 1, 1)
    histogram(diff_all)
    xlabel("FretMax - FretMin")
    
    subplot(2, 1, 2)
    histogram(diffNorm_all)
    xlabel("(FretMax - FretMin)/FretMin")
    
    %% Save outputs
    FretRange.diff_all = diff_all;
    FretRange.diffNorm_all = diffNorm_all;
    savefig(gcf, [OutputDest, 'FretRangeDistribution.fig'])
    saveas(gcf, [OutputDest, 'FretRangeDistribution.png'])
    
end