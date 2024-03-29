function fitSCHillFunction(popAvgData, SCAvgData, dataparms, parms, OutputDest)
    %% Make output directory
    mkdir([OutputDest, 'SingleCellHillFits/']); 

    %% Parameters
    R_norm = SCAvgData.ScMeanRnorm;
    R_norm_err = SCAvgData.ScStdRnorm;
    backConc = dataparms.backConc; %Background concentration of experiment
    concLevels = dataparms.concLevels; %Attractant concentrations measured
    Lplot = dataparms.Lplot;

    %% Fit Hill Function TO Each Cell
    %Equation for decreasing hill function. p(1) = log n, p(2) = log Kd, p(3) =
    %scaling factor
    %A*((1+((L-L0)/K)^n)^-1 - 1)
%     y = @(p, x) p(3) * (1./ (1 + ((x-backConc)./(exp(p(2)))).^exp(p(1))) - 1);
    y = @(p, x) (1./ (1 + ((x-backConc)./(exp(p(2)))).^exp(p(1))) - 1);

    if dataparms.Removal == 1
        y = @(p, x) p(3) * (1./ (1 + ((x)./(exp(p(2)))).^exp(p(1))));
    end

    p0 = dataparms.Hillp0; %initial guess for parameters
    
    for i = 1:size(R_norm, 1)
        DoseResp = R_norm(i, :);
%         sigma = R_norm_err(i, :);
        sigma = median(R_norm_err(i, :));
        inv_log_post = @(p) -(sum(-(1./(2*sigma.^2)) .* (DoseResp - y(p, concLevels)).^2)) - ...
            log(unifpdf(p(3), 0, 2));
        options = optimoptions('fminunc', 'MaxIterations', 10000, 'StepTolerance', 1e-07);
        p_opt(i, :) = fminunc(inv_log_post, p0, options);
    end
    
    %% Plots
    nPlots = 10; %Number of plots per figure
    for i = 1:size(R_norm, 1)        
        if mod(i, nPlots) == 1
            figure()
        end
        subplot(nPlots/2, 2, mod(i-1, nPlots) + 1)
        hold on
        plot(Lplot, y(p_opt(i, :), Lplot), 'r', 'Linewidth', 2)
        errorbar(concLevels, R_norm(i, :), R_norm_err(i, :), 'bo')
        set(gca, 'xscale', 'log')
        if ~dataparms.Removal
            ylim([-1.5, 0.2])
        else
            ylim([-0.5, 1.2])
        end
        
        xlim([(min(Lplot)), (max(Lplot))])
        
        if mod(i-1, nPlots) + 1 == nPlots | mod(i-1, nPlots) + 1 == nPlots -1
            xlabel('[L]')
        else
            set(gca, 'xtick', []) 
        end
        
        if mod(i-1, nPlots)+1 == nPlots | i == size(R_norm, 1) 
            set(gcf, 'Position', [100,25, 1200, 750])
            savefig(gcf, [OutputDest, 'SingleCellHillFits/', 'upto', num2str(i) ,'.fig'])
            saveas(gcf, [OutputDest, 'SingleCellHillFits/', 'upto', num2str(i) ,'.png'])
        end
    end
    
    %% Distributions of parameters
    %p_opt = (ncell, [log n, log K_1/2, A])
    Khalf = exp(p_opt(:,2));
    n = exp(p_opt(:, 1));
    A = p_opt(:, 3);
    
    %Filter out extreme fits
    Khalf(Khalf > median(Khalf) + std(Khalf)) = [];
    n(n > median(n) + std(n)) = [];
    A(A > median(A) + std(A)) = [];
    
    figure()
    subplot(3, 1, 1)
    histogram(Khalf)
    title("K_{1/2}, median = ", num2str(median(Khalf)));
    
    subplot(3, 1, 2)
    histogram(n)
    title("n, median = ", num2str(median(n)))
    
    subplot(3, 1, 3)
    histogram(A)
    title("A, median = ", num2str(median(A)))
    
    %% Parameter CDF
    cellRank = [1:length(Khalf)]./length(Khalf);
    Khalf_sorted = sort(Khalf);
    
    figure(); hold on
    plot(Khalf_sorted, cellRank, 'bo', 'linewidth', 2)
    set(gca, 'xscale', 'log')
    xlabel("K_{1/2}")
    ylabel("Cell Rank")
    title(['Median K_{1/2} = ', num2str(median(Khalf)), '\pm', num2str(std(Khalf))])
    savefig(gcf, [OutputDest, 'SingleCellHillFits/', 'RawCDF','.fig'])
    saveas(gcf, [OutputDest, 'SingleCellHillFits/', 'RawCDF','.png'])

end