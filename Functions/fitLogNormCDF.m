function CDFPlotData = fitLogNormCDF(SCAvgData, popAvgData, dataparms, parms, OutputDest)
    %% Collect points on the CDF using P(K<L) ~ P(R>0.5)
    concLevels = transpose(unique(dataparms.concLevels));
    concLevels = concLevels(:);
    
    %% Bootstrap single-cell average responses to get std for each point  
    %Calculate mean percent of cells that respond less than half-maximally
    singleCellRnorm = SCAvgData.ScMeanRnorm;
    NCells = sum(~isnan(singleCellRnorm)); %Number of cells for each concentration, ignores nans
    lessThanHalf = sum(singleCellRnorm < 0.5, 1);
%     if dataparms.Removal == 1
% %        singleCellRnorm = -singleCellRnorm; 
%        lessThanHalf = sum(singleCellRnorm < 0.5, 1); 
%     else 
%         lessThanHalf = sum(singleCellRnorm < -0.5, 1); 
%     end

    % Calculate using post-stimulus activity
    
    %Calculate confidence intervals and median with bootstrapping responses
    %in each cell
    nRuns = parms.nBootstraps;
%     scResponses = popAvgData.allRData_norm; %Responses normalized by prestimulus activity
%     scResponses = popAvgData.allRData; %Data not normalized by pre-stimulus
    scResponses = popAvgData.allPostStim;
    AllLessThanHalf = zeros(nRuns, size(lessThanHalf, 2));
    for i = 1:nRuns
        %First, resample each cells responses into a new dataset
        resamp = zeros(size(scResponses));
        
        %Pick random cells
        cellChoices = randi(size(scResponses, 3), size(scResponses, 3), 1);
        %Go through each chosen cell
        for jj = 1:size(scResponses, 3)
            c = cellChoices(jj);
            %Go through each stimulus level
            for s = 1:size(scResponses, 2)
                resamp(:, s, jj) = scResponses(randi(size(scResponses, 1), size(scResponses, 1),1), s, c);
            end
        end
        scA0 = popAvgData.a0(cellChoices);
        
        %Calculate median response for each cell
        med_resamp = squeeze(median(resamp, 1)); %(response, cell)
%         med_resamp = med_resamp ./ scA0; %Divide raw responses by a0. Don't use if calculating with fractional responses
%         med_resamp = (med_resamp - scA0) ./ scA0;
        med_resamp = (med_resamp) ./ scA0;

        %Calculate CDF points
%            if dataparms.Removal == 1
%                 AllLessThanHalf(i, :) = sum(med_resamp < 0.5, 2);
%            else
%                 AllLessThanHalf(i, :) = sum(med_resamp < -0.5, 2);
%            end
        AllLessThanHalf(i, :) = sum(med_resamp < 0.5, 2);
    end
%     CDFPoints = median(AllLessThanHalf, 1)./NCells;
    CDFPoints = lessThanHalf ./ NCells;
    CDFPointSamples = AllLessThanHalf ./ NCells;
    CDFPointErr = std(CDFPointSamples, 0, 1);
    CDFPointPos = abs(prctile(CDFPointSamples, 97.5) - CDFPoints);
    CDFPointNeg = abs(prctile(CDFPointSamples, 2.5) - CDFPoints);

    %% Bootstrap single-cell individual responses to get std for each point
    %Calculate mean percent of cells that respond less than half-maximally
%     IndividualRnorm = SCAvgData.allRData_norm; 
%     lessThanHalf = sum(IndividualRnorm < -0.5, 1); 
%     CDFPoints = lessThanHalf./size(IndividualRnorm, 1);
% 
%     %Calculate confidence intervals on each point with bootstrapping
%     nRuns = parms.nBootstraps;
%     AllLessThanHalf = zeros(nRuns, size(lessThanHalf, 2));
%     for i = 1:nRuns
%         %Get random samples from individual response data
%         samp = randi(size(IndividualRnorm, 1), size(IndividualRnorm, 1), size(IndividualRnorm, 2));
%         for j = 1:size(IndividualRnorm, 2)
%             ssRnormSamp(:, j) = IndividualRnorm(samp(:, j), j);
%         end
%         AllLessThanHalf(i, :) = sum(ssRnormSamp < -0.5, 1);
%     end
%     CDFPointSamples = AllLessThanHalf ./ size(IndividualRnorm, 1);
%     CDFPointErr = std(CDFPointSamples, 0, 1);
%     CDFPointPos = abs(prctile(CDFPointSamples, 97.5) - CDFPoints);
%     CDFPointNeg = abs(prctile(CDFPointSamples, 2.5) - CDFPoints);
   
    
    %% Fit Log-normal CDF to data
    sigma = CDFPointErr(:);
    CDFPoints = CDFPoints(:); %Reformat for fitting
    sigma(sigma == 0) = 0.0001; %Just in case you have a sample with 0 error
    p0 = dataparms.Normp0;
%     p0(1:2) = log(p0(1:2));
%     fitfun = @(b, L) b(3)*logncdf(L, (b(1)), ((b(2)))); %b(1) = meanKhalf, b(2) = variance, b(3) scaling
    fitfun = @(b, L) 1*logncdf(L, (b(1)), ((b(2)))); %b(1) = meanKhalf, b(2) = variance, b(3) scaling

    inv_log_post = @(p) -(sum(-(1./(2*sigma.^2)).*(CDFPoints - fitfun(p, (concLevels))).^2));
    p_opt2 = fminunc(inv_log_post, p0);
    
    % Log posterior with priors
    %Current priors are that:
        %n ~ N(n_ML, n_ML) (CV = n_ML)
        %K ~ N(K_ML, K_ML) (CV = 1)
        %A ~ U(0, 2)
    log_post = @(p) (sum(-(1./(2*sigma.^2)).*(CDFPoints - fitfun(p, (concLevels))).^2)) +...
        log(normpdf(p(1), p_opt2(1), abs(p_opt2(1)))) + log(normpdf(p(2), p_opt2(2), abs(p_opt2(2))));

    %% Evaluate posterior distribution to find error bars on parameters
    rnd = slicesample(p_opt2, parms.nSamples, 'logpdf', log_post, 'thin', parms.thining);
    
    
    %% Calculate percentiles
    %Convert samples to original units
    %Percentile for n
    K_low = prctile(rnd(:, 1), 2.5);
    K_high = prctile(rnd(:, 1), 97.5);
    Var_low = prctile(rnd(:, 2), 2.5);
    Var_high = prctile(rnd(:, 2), 97.5);
    
    %Find error bars
    K_err = (K_high - K_low) ./ 2;
    Var_err = (Var_high - Var_low)./ 2;
   
        %% Plot MCMC samples
    figure(); hold on

    subplot(3, 2, 1)
    plot(rnd(:, 1))
    ylabel("log K_{1/2}")
    subplot(3, 2, 2)
    hold on
    histogram(rnd(:, 1), 'Normalization', 'pdf')
    line([p_opt2(1), p_opt2(1)], ylim, 'LineWidth', 1, 'Color', 'r');
    camroll(-90)
    
    subplot(3, 2, 3)
    plot(rnd(:, 2))
    ylabel("log \sigma^2")
    subplot(3, 2, 4)
    histogram(rnd(:, 2), 'Normalization', 'pdf')
    line([p_opt2(2), p_opt2(2)], ylim, 'LineWidth', 1, 'Color', 'r');
    camroll(-90)

    
    %Save Figure in Output Destination
    savefig(gcf, [OutputDest, 'CDFMCMC.fig'])
    saveas(gcf, [OutputDest, 'CDFHillMCMC.png'])
    saveas(gcf, [OutputDest, 'CDFHillMCMC.svg'])


    
    %% Plot CDF Fit
%     p_opt2 = median(rnd, 1);
    
    %Extract plot range from parameters files
    LplotOld = dataparms.Lplot;
      
    if dataparms.backConc ~= 0
        Lplot = 10.^[log10(dataparms.backConc)-1:0.01:log10(max(LplotOld))];
    else
        Lplot = LplotOld;
    end
    
    %Generate individual curves for each of the MCMC Samples and calculate
    %95-CI
    curveSamples = zeros(size(rnd,1), length(Lplot));
    for i = 1:size(rnd, 1)
        curveSamples(i, :) = fitfun(rnd(i, :), Lplot);
    end
    curvePosErr = prctile(curveSamples, 97.5);
    curveNegErr = prctile(curveSamples, 2.5);
    
    %Plot figure
    figure('visible', parms.showPlot)
    hold on
        %Function to plot error of the curve as a fill using MCMC samples
        shadeLineError(Lplot, curvePosErr, curveNegErr, parms.lineColor)
    plot(Lplot, fitfun(p_opt2, Lplot), 'linewidth', 3, 'Color', parms.lineColor)
    errorbar(concLevels, CDFPoints, CDFPointNeg, CDFPointPos, 'o', 'MarkerSize', 5, 'Color', parms.lineColor)
    xlim([min(Lplot), max(Lplot)])
    set(gca, 'xscale', 'log')
    
    %Axis labels
    xlabel(dataparms.xlabels)
    ylabel("P(K_{1/2} < [L])")
    xlim([min(Lplot), max(Lplot)])
    ylim([0, 1.2])
    
    %Put data in title
    p_text = p_opt2;
%     title(['\langle K_{1/2} \rangle = ', num2str(p_text(1)), '\pm', num2str(K_err), '|', ...
%         '\sigma^2 = ', num2str(p_text(2)), '\pm', num2str(Var_err)])
    title(['log \langle K_{1/2} \rangle = ', num2str(p_text(1), 3), '[', num2str(K_low, 3), ',', num2str(K_high, 3),']', '|', ...
        '\sigma^2 = ', num2str(p_text(2), 3), '[', num2str(Var_low, 3), ',', num2str(Var_high, 3),']'])
    
    %Save Figure in Output Destination
    savefig(gcf, [OutputDest, 'CDFFit.fig'])
    saveas(gcf, [OutputDest, 'CDFFit.png'])
    saveas(gcf, [OutputDest, 'CDFFit.svg'])
    
    %% Save optimal parameters
    CDFPlotData.p_opt = p_opt2;
    CDFPlotData.Lplot = Lplot;
    CDFPlotData.concLevels = concLevels;
    CDFPlotData.CDFPoints = CDFPoints;
    CDFPlotData.K_err = K_err;
    CDFPlotData.Var_err = Var_err;
    CDFPlotData.CDFPointStd = CDFPointErr;
    CDFPlotData.CDFPointPosErr = CDFPointPos;
    CDFPlotData.CDFPointNegErr = CDFPointNeg;
    CDFPlotData.MCMCSamples = rnd;

end