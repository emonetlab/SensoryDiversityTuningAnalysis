%% Figure 4 supplement
% In this figure, we try to plot the mean meAsp K1/2 as a function of
% background L-asp. Ideally, and according to the model, this relationship
% should be linear. However, non-linearities in this relationship could
% explain some of the bizzarre features of figure 4, such as the inability
% to use the same ki values when Lasp or meAsp are in the background

close all; clear;
addpath('./Functions')

figure()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ====Extract K1/2 values===
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVFiles = ["./Outputs/MeAspDR_OneConcentrationSets/", ...
    "./Outputs/0_01LAsp_meAspDR/",...
     "./Outputs/1LAsp_meAspDR/",...
     "./Outputs/10LAspBackMeAspDose/",...
    ];

% metadata
backgroundLasp = [0, 0.01, 1, 10];

[CV, means, stds, CVerr, CVlowerbound, CVupperbound, meanerr,...
    meanlowbound, meanupperbound] = extract_cv_from_files(CVFiles);

plot(backgroundLasp, means, 'o-')
hold on
errorbar(backgroundLasp, means, means-meanlowbound, meanupperbound-means)
xlabel("[L-asp]")
ylabel("\langle meAsp K_{1/2}\rangle")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Khalf CDF w/ w/o Ser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data and declare plot parameters
% Files1 = {["./Outputs/0BackLAspDose/", "./Outputs/10serBackLAspDose/"]};
Files1 = {["./Outputs/0BackSerDose/", "./Outputs/10serBackLAspDose/"]};

% Color Palette To Color Each plot. Rows organized the same way
colrs = {["#4E4E4E", "#008080", "#008080", "#e02459",  "#e02459"]};

lnstyles = {["-", "--", "-", "--", "-"]};

CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));

plotInds = [1, 3];

xlims = [[10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3]];

xlabels = ["[L] \muM", "", "[meAsp] \muM", ""];

figure()
for ii = 1:length(Files1)
    FilesTemp = Files1{ii};
    colr = colrs{ii};
    xlimits = xlims(ii, :);
    lnstyl = lnstyles{ii};

    for jj = 1:length(FilesTemp)
        %Load Data
        load([convertStringsToChars(FilesTemp(jj)), 'plotData.mat'])
        Lplot = CDFPlotData.Lplot;

        %Generate sample curves for error bars
        rnd = CDFPlotData.MCMCSamples; %Dimensions: [mu, sigma, A (not used)]
        

        %Plot PDF on top
        subplot(4, 2, plotInds(ii, 1))
        hold on
        LplotPDF = log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot));
        y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
        plot(LplotPDF, y ./ max(y), 'color', colr(jj), 'Linewidth', 1.8,...
            'LineStyle', lnstyl(jj))
%         ylabel(['normalized', newline, 'pdf'])
        xlim([log(xlimits(1)), log(xlimits(2))])
        set(gca, 'XTick', [])
        if plotInds(ii, 2) /2 == round(plotInds(ii, 2)/2)
            set(gca, 'YTick', [])
        end
        if plotInds(ii, 2) /2 ~= round(plotInds(ii, 2)/2)
            ylabel('pdf')
        end
        if(jj == length(FilesTemp))
            set(gca, 'Position', get(gca, 'Position') + [0, -0.02, 0, -0.1])
        end
        

        %Plot CDF on bottom
        curveSamples = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples(j, :) = CDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr = prctile(curveSamples, 97.5);
        curveNegErr = prctile(curveSamples, 2.5);

        subplot(4, 2, plotInds(ii, 2))
        hold on
        shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(jj))
        errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr,...
            CDFPlotData.CDFPointPosErr, 'o', 'color', colr(jj));
        C(jj) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), ...
            'color', colr(jj), 'Linewidth', 1.8, 'LineStyle', lnstyl(jj));
        xlabel(xlabels(ii))
        xlim(xlimits)
        
        set(gca, 'xscale', 'log', 'XTick', [10^-2, 10^-1, 10^0, 10^1, 10^2])
        if plotInds(ii, 2) /2 ~= round(plotInds(ii, 2)/2)
            ylabel('P(K_{1/2} < [L])')
        else
            set(gca, 'YTick', [])
        end
        
    end
end

set(gcf, 'Position', [384,0,844,898])

%% Supplement part 2:
% Here we try to plot the population average dose-response curves for meAsp, and fit
% them using the MWC model to get an effective Ki and n.

Files = ["./Outputs/MeAspDR_OneConcentrationSets/plotData.mat",...
    "./Outputs/100MeAspBackMeAspDR/plotData.mat", ...
    "./Outputs/1meAspBackMeAspDR/plotData.mat",...
    "./Outputs/10MeAspBackMeAspDR/plotData.mat"];
bckg = [0, 100, 1, 10];
colr = ["b", "g", "c", "m", "r", "b"];

% Fitting function
% p = [n, e0, Ki]
f = @(L, L0, p) 1./(1 +((1+L./p(3))./(1+L0./p(3))).^p(1).*exp(p(2)));
% p0 = [40, 1, 60];
p0 = [12, 1, 18];

figure()
subplot(2, 1, 1)
x_store = zeros(length(Files), 5); 
y_store = zeros(length(Files), 5);

for ii = 1:length(Files)
    load(Files(ii))
    x = HillPlotData.concLevels;
    x_store(ii, :) = x;
    y = (HillPlotData.PopAvgDoseResp).*0.3;
    y_store(ii, :) = y;
%     y = (HillPlotData.PopAvgDoseResp);
%     yerr = (HillPlotData.PopAvgDoseStd);
end

xlim([10.^-1, 10.^3])
set(gca,'xscale', 'log')
xlabel("[MeAsp]")
ylabel("Kinase activity (a)")

% Fit values
inv_log_post = @(p) sum((y_store(1,:)-f(x_store(1,:), bckg(1), p)).^2) +...
    sum((y_store(2,:)-f(x_store(2,:), bckg(2), p)).^2) +...
    sum((y_store(3,:)-f(x_store(3,:), bckg(3), p)).^2) +...
    sum((y_store(4,:)-f(x_store(4,:), bckg(4), p)).^2);
p_opt = fminunc(inv_log_post, p0)

for ii = 1:length(Files)
    load(Files(ii))
    x = HillPlotData.concLevels;
    y = (HillPlotData.PopAvgDoseResp).*0.3;

    xfit = 10.^[-1:0.01:3];
    hold on
    plot(x, y, 'o', 'Color', colr(ii))
    plot(xfit, f(xfit, bckg(ii), p_opt), 'Color', colr(ii))
%     errorbar(x, y, yerr, 'o')
    hold off
end


% Try to predict the l-asp Khalf with meAsp present
CompFiles = ["./Outputs/10MeAsp10LAspBackMeAspDR/plotData.mat", ...
    "./Outputs/10LAsp_1meAspBackMeAspDR/plotData.mat",...
    "./Outputs/10LAsp_100meAsp_MeAspDR/plotData.mat",...
    "./Outputs/2LAsp_1MeAspBack_MeAspDR/plotData.mat",...
    "./Outputs/2LAsp_10MeAspBack_MeAspDR/plotData.mat"];
% p = [n, e0, Ki]
f_comp = @(L, L0, Lasp, p) 1./(1 +((1+L./p(3) + Lasp/0.85)./(1+L0./p(3) + Lasp/0.85)).^p(1).*exp(p(2)));
Lasp_bckg = [10, 10, 10, 2, 2];
meAsp_bckg = [10, 1, 100, 1, 10];

subplot(2, 1, 2)
for ii = 1:length(CompFiles)
    load(CompFiles(ii))
    x = HillPlotData.concLevels;
    y = (HillPlotData.PopAvgDoseResp).*0.3;

    xfit = 10.^[-2:0.01:3];
    hold on
    plot(x, y, 'o', 'Color', colr(ii))
    plot(xfit, f_comp(xfit,  meAsp_bckg(ii), Lasp_bckg(ii), p_opt), 'Color', colr(ii))
%     errorbar(x, y, yerr, 'o')
    hold off
end
xlabel("[MeAsp]")
ylabel("Kinase activity (a)")
set(gca, 'xscale', 'log')

%% Supplement part 2: L-asp responses
% Here we try to plot the population average dose-response curves for LAsp, and fit
% them using the MWC model to get an effective Ki and n.

Files = ["./Outputs/0BackLAspDose/plotData.mat",...
    "./Outputs/1LAsp_LAspDR/plotData.mat", ...
    "./Outputs/10LAsp_LAspDR/plotData.mat",...
    "./Outputs/0_1LAsp_LAspDR/plotData.mat"];
bckg = [0, 1, 10, 0.1];
colr = ["b", "g", "c", "m"];

% Fitting function
% p = [n, e0, Ki]
f = @(L, L0, p) 1./(1 +((1+L./p(3))./(1+L0./p(3))).^p(1).*exp(p(2)));
% p0 = [40, 1, 60];
p0 = [12, 1, 1];

figure()
subplot(2, 1, 1)
x_store = zeros(length(Files), 5); 
y_store = zeros(length(Files), 5);

for ii = 1:length(Files)
    load(Files(ii))
    x = HillPlotData.concLevels;
    x_store(ii, :) = x;
    y = (HillPlotData.PopAvgDoseResp).*0.3;
    y_store(ii, :) = y;
%     y = (HillPlotData.PopAvgDoseResp);
%     yerr = (HillPlotData.PopAvgDoseStd);
end

xlim([10.^-2, 10.^3])
set(gca,'xscale', 'log')
xlabel("[LAsp]")
ylabel("Kinase activity (a)")

% Fit values
inv_log_post = @(p) sum((y_store(1,:)-f(x_store(1,:), bckg(1), p)).^2) +...
    sum((y_store(2,:)-f(x_store(2,:), bckg(2), p)).^2) +...
    sum((y_store(3,:)-f(x_store(3,:), bckg(3), p)).^2) +...
    sum((y_store(4,:)-f(x_store(4,:), bckg(4), p)).^2);
p_opt = fminunc(inv_log_post, p0)

for ii = 1:length(Files)
    load(Files(ii))
    x = HillPlotData.concLevels;
    y = (HillPlotData.PopAvgDoseResp).*0.3;

    xfit = 10.^[-2:0.01:3];
    hold on
    plot(x, y, 'o', 'Color', colr(ii))
    plot(xfit, f(xfit, bckg(ii), p_opt), 'Color', colr(ii))
%     errorbar(x, y, yerr, 'o')
    hold off
end

% Try to predict the l-asp Khalf with meAsp present
CompFiles = ["./Outputs/100meAsp_LAspDR/plotData.mat", ...
    "./Outputs/100meAsp_0_1LAsp_LAspDR/plotData.mat",...
    "./Outputs/1LAsp_100meAsp_LAspDR/plotData.mat"];
% p = [n, e0, Ki]
f_comp = @(L, L0, meAsp, p) 1./(1 +((1+L./p(3) + meAsp/17.2)./(1+L0./p(3) + meAsp/17.2)).^p(1).*exp(p(2)));
Lasp_bckg = [0, 0.1, 1];
meAsp_bckg = [100, 100, 100];

subplot(2, 1, 2)
for ii = 1:length(CompFiles)
    load(CompFiles(ii))
    x = HillPlotData.concLevels;
    y = (HillPlotData.PopAvgDoseResp).*0.3;

    xfit = 10.^[-2:0.01:3];
    hold on
    plot(x, y, 'o', 'Color', colr(ii))
    plot(xfit, f_comp(xfit, Lasp_bckg(ii), meAsp_bckg(ii), p_opt), 'Color', colr(ii))
%     errorbar(x, y, yerr, 'o')
    hold off
end
xlabel("[LAsp]")
ylabel("Kinase activity (a)")
set(gca, 'xscale', 'log')
