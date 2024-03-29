close all; clear;
addpath('./Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Organizing Data for fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Piece 1: plot the CV curves. 1mMGluBackMeAspDose MeAspDR_OneConcentrationSets
CVFiles = ["./Outputs/MeAspDR_OneConcentrationSets/", ...
    "./Outputs/10MeAspBackMeAspDR/",...
    "./Outputs/100MeAspBackMeAspDR/",...
    "./Outputs/0_3MeAspBackMeAspDR/",...
    "./Outputs/0_1meAspBackMeAspDR/",...
    "./Outputs/1meAspBackMeAspDR/",...
    "./Outputs/0_01meAsp_meAspDR/",...
    "./Outputs/10MeAsp10LAspBackMeAspDR/",...
    "./Outputs/10LAsp_1meAspBackMeAspDR/",...
    "./Outputs/10LAsp_100meAsp_MeAspDR/",...
    "./Outputs/2LAsp_1MeAspBack_MeAspDR/",...
    "./Outputs/2LAsp_10MeAspBack_MeAspDR/",...
    "./Outputs/0BackLAspDose/", ... % Begin L_asp dose-responses
    "./Outputs/1LAsp_LAspDR/", ...
    "./Outputs/0_01LAsp_LAspDR/", ...
    "./Outputs/0_5LAsp_LAspDR/",...
    "./Outputs/0_1LAsp_LAspDR/", ...
    "./Outputs/0_05LAsp_LAspDR/",...
    "./Outputs/10LAsp_LAspDR/",...
    "./Outputs/1LAsp_100meAsp_LAspDR/",...
    "./Outputs/10LAsp_100meAsp_LAspDR_V2/",...
    "./Outputs/100MeAspBackLAspDose/",...
    "./Outputs/100meAsp_0_1LAsp_LAspDR/", ...
    "./Outputs/1MaltosebackMeAspDose/",...
    "./Outputs/10MeAsp10LAspBackMeAspDR/"];

% Gave 0 background small value for plot, but fit with 0
bckg = [0, 10, 100, 0.3, 0.1, 1, 0.01, 10, 1, 100, 1, 10,...
    0, 0, 0, 0, 0, 0, 0, 100, 100, 100, 100, 0];
bckg_plot = bckg; bckg_plot(bckg==0) = 10^-4;  
LAspBckg = [0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 2, 2,...
    0, 1, 0.01, 0.5, 0.1, 0.05, 10, 1, 10, 0, 0.1, 0];
bckg_plot_0 = LAspBckg; bckg_plot_0(bckg_plot_0==0) = 10^-4;  
ligandchoice = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1];

% Extract CV from each file, and store in order in vectors
% [CV, means, stds, err] = extract_cv_from_files(CVFiles);
[CV, means, stds, CVerr, CVlowerbound, CVupperbound, meanerr,...
    meanlowbound, meanupperbound] = extract_cv_from_files(CVFiles);


%Calculate FCKhalf curve as a function of background
colr = ["#4E4E4E"; "#798086"; "#50B2C0"];
mrkrSize = 7;
% 
L0Range = 10.^[-4:0.1:3];

% Calculate fitting parameters for meAsp khalf data
backgrounds = bckg(LAspBckg==0 & ligandchoice==1); 
backgrounds2 = bckg(LAspBckg==2 & ligandchoice==1); 
backgrounds10 = bckg(LAspBckg==10 & ligandchoice==1);
y_fit = CV(LAspBckg==0 & ligandchoice==1); 
y_fit2 = CV(LAspBckg==2 & ligandchoice==1); 
y_fit10 = CV(LAspBckg==10& ligandchoice==1);
err_fit = CVerr(LAspBckg==0& ligandchoice==1); 
err_fit2 = CVerr(LAspBckg==2& ligandchoice==1); 
err_fit10 = CVerr(LAspBckg==10& ligandchoice==1);

% Calculate fitting parameters for Lasp khalf data
backgrounds_0 = LAspBckg(bckg==0& ligandchoice==0);
y_fit_0 = CV(bckg==0& ligandchoice==0);
err_fit_0 = CVerr(bckg==0& ligandchoice==0);
backgrounds2_0 = LAspBckg(bckg==100& ligandchoice==0);
y_fit2_0 = CV(bckg==100& ligandchoice==0);
err_fit2_0 = CVerr(bckg==100& ligandchoice==0);

% means(1)
% stds(1)
% means(13)
% stds(13)
CV(end)
CVerr(end)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select model to fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Model with full negative cooperativity
% CVfunc = @negative_cooperativity_model;
% % p0 = [Ki1, Ki2, e11, e12, e22, mu, sigma]
% p0 = [5, 18, 6, 3, 5, -1, 2];
% inv_log_post = @(p) sum((y_fit - CVfunc(p, 0, backgrounds, 1)).^2./err_fit.^2) +... 
%     sum((y_fit2 - CVfunc(p, 2, backgrounds2, 1)).^2./err_fit2.^2) + ...
%     sum((y_fit10 - CVfunc(p, 10, backgrounds10, 1)).^2./err_fit10.^2) + ...
%     sum((y_fit_0 - CVfunc(p, backgrounds_0, 0, 0)).^2./err_fit_0.^2) + ...
%     sum((y_fit2_0 - CVfunc(p, backgrounds2_0, 100, 0)).^2./err_fit2_0.^2) + ... % Begin priors
%     unifpdf(p(1), 0, 100) + unifpdf(p(2), 0, 100) + unifpdf(p(3), 0, 20) +...
%     unifpdf(p(4), 0, 20) + unifpdf(p(5), 0, 12) +...
%     unifpdf(p(6), -10, 10) + unifpdf(p(7), 0, 5);
% 
% log_post = @(p) -(sum((y_fit - CVfunc(p, 0, backgrounds, 1)).^2./err_fit.^2) +... 
%     sum((y_fit2 - CVfunc(p, 2, backgrounds2, 1)).^2./err_fit2.^2) + ...
%     sum((y_fit10 - CVfunc(p, 10, backgrounds10, 1)).^2./err_fit10.^2) + ...
%     sum((y_fit_0 - CVfunc(p, backgrounds_0, 0, 0)).^2./err_fit_0.^2) + ...
%     sum((y_fit2_0 - CVfunc(p, backgrounds2_0, 100, 0)).^2./err_fit2_0.^2))+ ... % Begin priors
%     log(unifpdf(p(1), 0, 100)) + log(unifpdf(p(2), 0, 100)) + log(unifpdf(p(3), 0, 20)) +...
%     log(unifpdf(p(4), 0, 20))+ log(unifpdf(p(5), 0, 12)) +...
%     log(unifpdf(p(6), -10, 10)) + log(unifpdf(p(7), 0, 5));

% Model with two binding sites
% CVfunc = @two_binding_site_model;
% % p0 = [Ki1, Ki2, Ki11, mu, sigma]
% p0 = [0.8, 18, 1, -1, .1];
% inv_log_post = @(p) sum((y_fit - CVfunc(p, 0, backgrounds, 1)).^2./err_fit.^2) +... 
%     sum((y_fit2 - CVfunc(p, 2, backgrounds2, 1)).^2./err_fit2.^2) + ...
%     sum((y_fit10 - CVfunc(p, 10, backgrounds10, 1)).^2./err_fit10.^2) + ...
%     sum((y_fit_0 - CVfunc(p, backgrounds_0, 0, 0)).^2./err_fit_0.^2) + ...
%     sum((y_fit2_0 - CVfunc(p, backgrounds2_0, 100, 0)).^2./err_fit2_0.^2) + ... % Begin priors
%     log(normpdf(p(1), 0.8, 0.5)) + log(normpdf(p(2), 18, 1)) + ...
%     log(unifpdf(p(3), 0, 50)) + log(unifpdf(p(4), -10, 10)) +...
%     log(unifpdf(p(5), 0, 5));
% 
% log_post = @(p) -(sum((y_fit - CVfunc(p, 0, backgrounds, 1)).^2./err_fit.^2) +... 
%     sum((y_fit2 - CVfunc(p, 2, backgrounds2, 1)).^2./err_fit2.^2) + ...
%     sum((y_fit10 - CVfunc(p, 10, backgrounds10, 1)).^2./err_fit10.^2) + ...
%     sum((y_fit_0 - CVfunc(p, backgrounds_0, 0, 0)).^2./err_fit_0.^2) + ...
%     sum((y_fit2_0 - CVfunc(p, backgrounds2_0, 100, 0)).^2./err_fit2_0.^2))+ ... % Begin priors
%     log(unifpdf(p(1), 0, 100)) + log(unifpdf(p(2), 0, 100)) + log(unifpdf(p(3), 0, 100)) +...
%     log(unifpdf(p(4), -10, 10)) + log(unifpdf(p(5), 0, 5));
% 
% % Model with one binding site
% CVfunc = @one_binding_site_model;
% % p0 = [Ki1, Ki2, mu, sigma]
% p0 = [0.1, 5, -1, .1];
% inv_log_post = @(p) sum((y_fit - CVfunc(p, 0, backgrounds, 1)).^2./err_fit.^2) +... 
%     sum((y_fit2 - CVfunc(p, 2, backgrounds2, 1)).^2./err_fit2.^2) + ...
%     sum((y_fit10 - CVfunc(p, 10, backgrounds10, 1)).^2./err_fit10.^2) + ...
%     sum((y_fit_0 - CVfunc(p, backgrounds_0, 0, 0)).^2./err_fit_0.^2) + ...
%     sum((y_fit2_0 - CVfunc(p, backgrounds2_0, 100, 0)).^2./err_fit2_0.^2) + ... % Begin priors
%     log(unifpdf(p(1), 0, 100)) + log(unifpdf(p(2), 0, 100)) +...
%     log(unifpdf(p(3), -10, 10)) + log(unifpdf(p(4), 0, 5));
% 
% log_post = @(p) -(sum((y_fit - CVfunc(p, 0, backgrounds, 1)).^2./err_fit.^2) +... 
%     sum((y_fit2 - CVfunc(p, 2, backgrounds2, 1)).^2./err_fit2.^2) + ...
%     sum((y_fit10 - CVfunc(p, 10, backgrounds10, 1)).^2./err_fit10.^2) + ...
%     sum((y_fit_0 - CVfunc(p, backgrounds_0, 0, 0)).^2./err_fit_0.^2) + ...
%     sum((y_fit2_0 - CVfunc(p, backgrounds2_0, 100, 0)).^2./err_fit2_0.^2))+ ... % Begin priors
%     log(unifpdf(p(1), 0, 100)) + log(unifpdf(p(2), 0, 100)) +...
%     log(unifpdf(p(3), -10, 10)) + log(unifpdf(p(4), 0, 5));

% Deterministic Model with one binding site
CVfunc = @one_binding_site_model_deterministic;
% p0 = [Ki1, Ki2, K0_meAsp, sigma0_meAsp, K0_Lasp, sigma0_Lasp]
p0 = [0.85, 17.2, means(1), stds(1), means(13), stds(13)];
inv_log_post = @(p) (sum((y_fit - CVfunc(p, 0, backgrounds, 1)).^2) +... 
    1*sum((y_fit2 - CVfunc(p, 2, backgrounds2, 1)).^2) + ...
    1*sum((y_fit10 - CVfunc(p, 10, backgrounds10, 1)).^2) + ...
        1*(sum((y_fit_0 - CVfunc(p, backgrounds_0, 0, 0)).^2) + ...
    sum((y_fit2_0 - CVfunc(p, backgrounds2_0, 100, 0)).^2))- ... 
    1*(log(normpdf(p(1), 0.85, 0.1)) + log(normpdf(p(2), 18,1)) +...
    log(normpdf(p(3), 1.25, 0.2)) + log(unifpdf(p(4), 0, 10)) + ...
    log(normpdf(p(5), 0.0565, 0.005)) + log(unifpdf(p(6), 0, 10))))

log_post = @(p) -(sum((y_fit - CVfunc(p, 0, backgrounds, 1)).^2) +... 
    1*sum((y_fit2 - CVfunc(p, 2, backgrounds2, 1)).^2) + ...
    1*sum((y_fit10 - CVfunc(p, 10, backgrounds10, 1)).^2) + ...
        1*(sum((y_fit_0 - CVfunc(p, backgrounds_0, 0, 0)).^2) + ...
    sum((y_fit2_0 - CVfunc(p, backgrounds2_0, 100, 0)).^2))- ... 
    1*(log(normpdf(p(1), 0.85, 0.1)) + log(normpdf(p(2), 18,1)) +...
    log(normpdf(p(3), 1.25, 0.2)) + log(unifpdf(p(4), 0, 10)) + ...
    log(normpdf(p(5), 0.0565, 0.005)) + log(unifpdf(p(6), 0, 10))))



options = optimoptions("fminunc",StepTolerance=1e-09);
p_opt_temp = fminunc(inv_log_post, p0);

disp('Optimal CV Curve Parameters:')
p0
p_opt = p_opt_temp
% p_opt = p0


%% mcmc sampling
rnd = slicesample(p_opt_temp, 2000, "logpdf", log_post, "thin", 10);
p_samp = median(rnd);
p_samp
std(rnd)
figure()
hold on
for ii = 1:length(p_samp)
    subplot(length(p_samp), 1, ii)
    plot(rnd(:, ii))
end
% p0 = [Ki1, Ki2, K0_meAsp, sigma0_meAsp, K0_Lasp, sigma0_Lasp]
p_opt = median(rnd);
%% ===Plotting===
xlimits = [10.^-4.5, 10.^3];
figure()
subplot(4, 2, [1, 3])
CV0 = CVfunc(p_opt, 0, L0Range, 1);
CV10 = CVfunc(p_opt, 10, L0Range, 1);
CV2 = CVfunc(p_opt, 2, L0Range, 1);
% CV0 = CVfitfunc(p0, L0Range);
% CV10 = CVfitfunc_competitor(p0, L0Range, 10);
% CV2 = CVfitfunc_competitor(p0, L0Range, 2);

hold on
C(1) = plot(L0Range, CV0, 'Linewidth', 2, 'Color', colr(1,:));
C(2) = plot(L0Range, CV10, 'Linewidth', 2, 'Color', colr(3,:));
C(3) = plot(L0Range, CV2, 'Linewidth', 2, 'Color', colr(2,:));
errorbar(bckg_plot(LAspBckg==0 & ligandchoice==1), CV(LAspBckg==0 & ligandchoice==1),  CV(LAspBckg==0 & ligandchoice==1)-CVlowerbound(LAspBckg==0 & ligandchoice==1),...
    CVupperbound(LAspBckg==0& ligandchoice==1)-CV(LAspBckg==0& ligandchoice==1), 'o', 'MarkerFaceColor', colr(1,:), 'MarkerEdgeColor', colr(1,:),...
    'Color', colr(1,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
errorbar(bckg_plot(LAspBckg==10& ligandchoice==1), CV(LAspBckg==10& ligandchoice==1), CV(LAspBckg==10& ligandchoice==1)-CVlowerbound(LAspBckg==10& ligandchoice==1),...
    CVupperbound(LAspBckg==10& ligandchoice==1)-CV(LAspBckg==10& ligandchoice==1), 'o', 'MarkerFaceColor', colr(3,:), 'MarkerEdgeColor', colr(3,:),...
    'Color', colr(3,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
errorbar(bckg_plot(LAspBckg==2& ligandchoice==1), CV(LAspBckg==2& ligandchoice==1), CV(LAspBckg==2& ligandchoice==1)-CVlowerbound(LAspBckg==2& ligandchoice==1),...
    CVupperbound(LAspBckg==2& ligandchoice==1)-CV(LAspBckg==2& ligandchoice==1), 'o', 'MarkerFaceColor', colr(2,:), 'MarkerEdgeColor', colr(2,:),...
    'Color', colr(2,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
ylabel("K_{1/2} CV")
xlabel("Background [meAsp] \muM")
xlim(xlimits)
% ylim([0,.75])
xticks([10.^-4, 10.^-3, 10^-2, 10^-1, 10^-0, 10^1, 10^2])
xticklabels(["0", "10^{-3}", "10^{-2}", "10^{-1}", "10^0", "10^1", "10^2"])
set(gca, 'xscale', 'log', 'XMinorTick', 'off')
hold off

subplot(4, 2, [5, 7])
CV0 = CVfunc(p_opt, L0Range, 0, 0);
CV100 = CVfunc(p_opt, L0Range, 100, 0);
hold on
C(1) = plot(L0Range, CV0, 'Linewidth', 2, 'Color', colr(1,:));
C(2) = plot(L0Range, CV100, 'Linewidth', 2, 'Color', colr(2,:));
errorbar(bckg_plot_0(bckg==0& ligandchoice==0), CV(bckg==0& ligandchoice==0), CV(bckg==0& ligandchoice==0)-CVlowerbound(bckg==0& ligandchoice==0),...
    CVupperbound(bckg==0& ligandchoice==0)-CV(bckg==0& ligandchoice==0), 'o', 'MarkerFaceColor', colr(1,:), 'MarkerEdgeColor',  colr(1,:),...
    'Color', colr(1,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
errorbar(bckg_plot_0(bckg==100& ligandchoice==0), CV(bckg==100& ligandchoice==0), CV(bckg==100& ligandchoice==0)-CVlowerbound(bckg==100& ligandchoice==0),...
    CVupperbound(bckg==100& ligandchoice==0)-CV(bckg==100& ligandchoice==0), 'o', 'MarkerFaceColor', colr(2,:), 'MarkerEdgeColor',  colr(2,:),...
    'Color', colr(2,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
ylabel("K_{1/2} CV")
xlabel("Background [L-Asp] \muM")
xlim(xlimits)
ylim([0, 1.7])
xticks([10.^-4, 10.^-3, 10^-2, 10^-1, 10^-0, 10^1, 10^2])
xticklabels(["0", "10^{-3}", "10^{-2}", "10^{-1}", "10^0", "10^1", "10^2"])
set(gca, 'xscale', 'log', 'XMinorTick', 'off')


%% Figure aesthetics
% set(gcf, 'Position', [151.4,97,896.6,665])
set(gcf, 'Position', [151.4,1,896.6,750])

saveas(gcf, './Outputs/Figure4_V5_part1_negcoop.svg')
saveas(gcf, './Outputs/Figure4_V5_part1_negcoop.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting meAsp khalf distribution with Lasp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Files are organized by row:
    %Row 1: top left
    %Row 2: bottom left
    %Row 3: top right
    %Row 4: bottom right
Files1 = {["./Outputs/MeAspDR_OneConcentrationSets/", ...
    "./Outputs/10MeAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/", ...
    "./Outputs/100MeAspBackMeAspDR/"], ...
    ["./Outputs/0BackLAspDose/",...
    "./Outputs/1LAsp_LAspDR/", ...
    "./Outputs/10LAsp_LAspDR/",...
    "./Outputs/1LAsp_100meAsp_LAspDR/",...
    "./Outputs/10LAsp_100meAsp_LAspDR_V2/"]};
% Files1 = {["./Outputs/0BackLAspDose/", ...
%     "./Outputs/10MeAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/"]};



% Color Palette To Color Each plot. Rows organized the same way
colrs = {["#4E4E4E", "#058ED9", "#e02459", "#058ED9", "#058ED9"], ...
    ["#4E4E4E", "#058ED9", "#e02459", "#058ED9", "#058ED9"]};

lnstyles = {["-", "-", "-", "--", "--", "--"], ["-", "-", "-", "--", "--", "--"]};

showMarkers = {[1, 1, 1, 0, 1], [1, 1, 1, 1, 1]};

backgrColors = {[1 1 1], [1 1 1]};

CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));


% Plot
%Establish grid: Should be 4x2

%Position of the plots organized in the same order as Files1
plotInds = [2, 4;
            6, 8];

xlims = [[10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3];...
       [10^-3, 10^3]];

xlabels = ["[meAsp] \muM", "", "[meAsp] \muM", ""];

CV = [];
means = [];
stds = [];
err = [];
% bckg = zeros(length(Files1));

figure()

% Plot PDF on top
for ii = 1:length(Files1)
    FilesTemp = Files1{ii};
    colr = colrs{ii};
    xlimits = xlims(ii, :);
    showMarker = showMarkers{ii};
    lnstyle = lnstyles{ii};

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
        plot(LplotPDF, y ./ max(y), 'color', colr(jj), 'Linewidth', 1.8, ...
           'LineStyle', lnstyle(jj))
%         ylabel(['normalized', newline, 'pdf'])
        xlim([log(xlimits(1)), log(xlimits(2))])
        set(gca, 'XTick', [], 'color', backgrColors{ii})
        ylabel('pdf')
%         if(jj == length(FilesTemp))
%             set(gca, 'Position', get(gca, 'Position') + [0, -0.02, 0, -0.1])
%         end

        

        %Plot CDF on bottom
        curveSamples = zeros(size(rnd,1), length(Lplot));
        for j = 1:size(rnd, 1)
            curveSamples(j, :) = CDFFunc(rnd(j, :), Lplot);
        end
        curvePosErr = prctile(curveSamples, 97.5);
        curveNegErr = prctile(curveSamples, 2.5);

        subplot(4, 2, plotInds(ii, 2))
        hold on
        if showMarker(jj)
            shadeLineError(CDFPlotData.Lplot, curvePosErr, curveNegErr, colr(jj))
            errorbar(CDFPlotData.concLevels, CDFPlotData.CDFPoints, CDFPlotData.CDFPointNegErr,...
                CDFPlotData.CDFPointPosErr, 'o', 'color', colr(jj));
        end
        C(jj) = plot(CDFPlotData.Lplot, CDFFunc(CDFPlotData.p_opt, CDFPlotData.Lplot), ...
            'color', colr(jj), 'Linewidth', 1.8, 'LineStyle', lnstyle(jj));
        xlabel(xlabels(ii))
        xlim(xlimits)

        set(gca, 'xscale', 'log', 'XTick', [10^-2, 10^-1, 10^0, 10^1, 10^2],...
            'color', backgrColors{ii})
        ylabel('P(K_{1/2} < [L])')
        
        CDFPlotData.p_opt;

        %Calculate CV
        mu_meAsp = CDFPlotData.p_opt(1);
        sigma = CDFPlotData.p_opt(2);

        m_all = exp(rnd(:, 1) + rnd(:,2).^2./2);
        v_all = exp(2*rnd(:, 1) + rnd(:, 2).^2).*(exp(rnd(:, 2).^2) - 1);

    end
end

set(gcf, 'Position', [384,2,844,898])

%% Plot predicted Khalf distributions
% Plot predicted khalf distributions using fitted meAsp model parameters
% p0 = [Ki_Lasp, Ki_meAsp, K0_meAsp, sigma0_meAsp, K0_Lasp, sigma0_Lasp]
params = p_opt;
CVfitfunc_Lasp = @(p, meAsp, Lasp) p(6) ./ (p(5) + Lasp./(1 + meAsp./p(2) + Lasp./p(1)));
Meanfunc_Lasp = @(p, meAsp, Lasp) p(5)*(1+meAsp./p(2) + Lasp./p(1)) + Lasp;
CVfitfunc_meAsp = @(p, meAsp, Lasp) p(4) ./ (p(3) + meAsp./(1 + meAsp./p(2) + Lasp./p(1)));
Meanfunc_meAsp = @(p, meAsp, Lasp) p(4)*(1+meAsp./p(2) + Lasp./p(1)) + meAsp;

% Get lognormal parameters for meAsp from mean and CV
m_meAsp = @(params, meAsp, Lasp) Meanfunc_meAsp(params, meAsp, Lasp);
v_meAsp = @(params, meAsp, Lasp) (m_meAsp(params, meAsp, Lasp).*CVfitfunc_meAsp(params, meAsp, Lasp)).^2;
mu_meAsp = @(params, meAsp, Lasp) log(m_meAsp(params, meAsp, Lasp).^2./(sqrt(v_meAsp(params, meAsp, Lasp)+m_meAsp(params, meAsp, Lasp).^2)));
sig_meAsp = @(params, meAsp, Lasp) sqrt(log(v_meAsp(params, meAsp, Lasp)./m_meAsp(params, meAsp, Lasp).^2 + 1));


Lrange = 10.^[-3:0.01:3];
subplot(4, 2, plotInds(1, 2))
hold on
meAsp = 0; Lasp = 0;
plot(Lrange, logncdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp)),  'k--', 'linewidth', 2)
meAsp = 10; Lasp = 0;
plot(Lrange, logncdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp)),  'k--', 'linewidth', 2)
meAsp = 10; Lasp = 10;
plot(Lrange, logncdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp)), 'k--', 'linewidth', 2)
meAsp = 100; Lasp = 10;
plot(Lrange, logncdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp)), 'k--', 'linewidth', 2)
set(gca, 'xscale', 'log', 'XTick', [10^-2, 10^-1, 10^0, 10^1, 10^2])
xlim([10.^-3, 10.^3])

subplot(4, 2, plotInds(1, 1))
hold on
Lrange = -4:0.01:4;
meAsp = 0; Lasp = 0;
plot(Lrange, normpdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp))./...
    max(normpdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp))),...
    'k--', 'linewidth', 2)

Lrange = 2:0.01:6;
meAsp = 10; Lasp = 0;
plot(Lrange, normpdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp))./...
    max(normpdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp))),...
    'k--', 'linewidth', 2)
Lrange = 2:0.01:6;
meAsp = 10; Lasp = 10;
plot(Lrange, normpdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp))./...
    max(normpdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp))),...
    'k--', 'linewidth', 2)

Lrange = 3:0.01:7;
meAsp = 100; Lasp = 10;
plot(Lrange, normpdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp))./...
    max(normpdf(Lrange, mu_meAsp(params, meAsp, Lasp), sig_meAsp(params, meAsp, Lasp))),...
    'k--', 'linewidth', 2)
set(gca, 'Position', get(gca, 'Position') + [0, -0.02, 0, -0.1])

%% Plot a grid of measured and predicted Khalf distributions
Files1 = ["./Outputs/MeAspDR_OneConcentrationSets/", ...
    "./Outputs/0_01meAsp_meAspDR/",...
    "./Outputs/0_1meAspBackMeAspDR/",...
    "./Outputs/0_3MeAspBackMeAspDR/",...
    "./Outputs/1meAspBackMeAspDR/",...
    "./Outputs/10MeAspBackMeAspDR/",...
    "./Outputs/100MeAspBackMeAspDR/",...
    "./Outputs/2LAsp_1MeAspBack_MeAspDR/",...
    "./Outputs/2LAsp_10MeAspBack_MeAspDR/",...
    "./Outputs/10LAsp_1meAspBackMeAspDR/",...
    "./Outputs/10MeAsp10LAspBackMeAspDR/",...
    "./Outputs/10LAsp_100meAsp_MeAspDR/"];
meAsp_list1 = [0, 0.01, 0.1, 0.3, 1, 10, 100, 1, 10, 1, 10, 100];
Lasp_list1 = [0, 0, 0, 0, 0, 0, 0, 2, 2, 10, 10, 10];


Files2 = [ "./Outputs/0BackLAspDose/", ... % Begin L_asp dose-responses
    "./Outputs/0_01LAsp_LAspDR/", ...
    "./Outputs/0_05LAsp_LAspDR/",...
    "./Outputs/0_1LAsp_LAspDR/", ...
    "./Outputs/0_5LAsp_LAspDR/",...
    "./Outputs/1LAsp_LAspDR/", ...
    "./Outputs/10LAsp_LAspDR/",...
    "./Outputs/100MeAspBackLAspDose/",...
    "./Outputs/100meAsp_0_1LAsp_LAspDR/",...
    "./Outputs/1LAsp_100meAsp_LAspDR/",...
    "./Outputs/10LAsp_100meAsp_LAspDR_V2/"];

meAsp_list2 = [0,0,0,0,0,0,0,100, 100, 100, 100];
Lasp_list2 = [0, 0.01, 0.05, 0.1, 0.5, 1, 10, 0, 0.1, 1, 10];

CDFFunc = @(p, L) logncdf(L, p(1), p(2));
PDFFunc = @(p, L) normpdf(L, p(1), p(2));

% Equations for predicted Khalf distribution
params = p_opt;
CVfitfunc_Lasp = @(p, meAsp, Lasp) p(6) ./ (p(5) + Lasp./(1 + meAsp./p(2) + Lasp./p(1)));
Meanfunc_Lasp = @(p, meAsp, Lasp) p(5)*(1+meAsp./p(2) + Lasp./p(1)) + Lasp;
CVfitfunc_meAsp = @(p, meAsp, Lasp) p(4) ./ (p(3) + meAsp./(1 + meAsp./p(2) + Lasp./p(1)));
Meanfunc_meAsp = @(p, meAsp, Lasp) p(4)*(1+meAsp./p(2) + Lasp./p(1)) + meAsp;
% Get lognormal parameters for meAsp from mean and CV
m_meAsp = @(params, meAsp, Lasp) Meanfunc_meAsp(params, meAsp, Lasp);
v_meAsp = @(params, meAsp, Lasp) (m_meAsp(params, meAsp, Lasp).*CVfitfunc_meAsp(params, meAsp, Lasp)).^2;
mu_meAsp = @(params, meAsp, Lasp) log(m_meAsp(params, meAsp, Lasp).^2./(sqrt(v_meAsp(params, meAsp, Lasp)+m_meAsp(params, meAsp, Lasp).^2)));
sig_meAsp = @(params, meAsp, Lasp) sqrt(log(v_meAsp(params, meAsp, Lasp)./m_meAsp(params, meAsp, Lasp).^2 + 1));

m_Lasp = @(params, meAsp, Lasp) Meanfunc_Lasp(params, meAsp, Lasp);
v_LAsp = @(params, meAsp, Lasp) (m_Lasp(params, meAsp, Lasp).*CVfitfunc_Lasp(params, meAsp, Lasp)).^2;
mu_Lasp = @(params, meAsp, Lasp) log(m_Lasp(params, meAsp, Lasp).^2./(sqrt(v_LAsp(params, meAsp, Lasp)+m_Lasp(params, meAsp, Lasp).^2)));
sig_Lasp = @(params, meAsp, Lasp) sqrt(log(v_LAsp(params, meAsp, Lasp)./m_Lasp(params, meAsp, Lasp).^2 + 1));



CV = [];
means = [];
stds = [];
err = [];
% bckg = zeros(length(Files1));



figure()
nrows = round(length(Files1)./2);
for jj = 1:length(Files1)
    %Load Data
    load([convertStringsToChars(Files1(jj)), 'plotData.mat'])

    Lplot = CDFPlotData.Lplot;

    %Generate sample curves for error bars
    rnd = CDFPlotData.MCMCSamples; %Dimensions: [mu, sigma, A (not used)]
    

    %Get error for PDF plot
        LplotPDF = log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot));
     curveSamples = zeros(size(rnd,1), length(LplotPDF));
    for j = 1:size(rnd, 1)
        curveSamples(j, :) = PDFFunc(rnd(j, :), LplotPDF);
    end
    curvePosErr = prctile(curveSamples, 97.5);
    curveNegErr = prctile(curveSamples, 2.5);

    subplot(nrows, 2, jj)
    hold on
    LplotPDF = log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot));
    y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
    shadeLineError(LplotPDF, curvePosErr, curveNegErr, [0.2, 0.2, 0.2])
    plot(LplotPDF, y, 'k-', 'Linewidth', 1)
%         ylabel(['normalized', newline, 'pdf'])
    xlim([log(xlimits(1)), log(xlimits(2))])
    ylabel('pdf')
%         if(jj == length(FilesTemp))
%             set(gca, 'Position', get(gca, 'Position') + [0, -0.02, 0, -0.1])
%         end
    meAsp_temp = meAsp_list1(jj); 
    Lasp_temp = (Lasp_list1(jj));
%     plot(LplotPDF, normpdf(LplotPDF, mu_meAsp(params, meAsp_temp, Lasp_temp), sig_meAsp(params, meAsp_temp, Lasp_temp))./...
%         max(normpdf(LplotPDF, mu_meAsp(params, meAsp_temp, Lasp_temp), sig_meAsp(params, meAsp_temp, Lasp_temp))),...
%         'r--', 'linewidth', 2)
    plot(LplotPDF, normpdf(LplotPDF, mu_meAsp(params, meAsp_temp, Lasp_temp), sig_meAsp(params, meAsp_temp, Lasp_temp)),...
        'r--', 'linewidth', 1)
end


set(gcf, 'Position', [384,-9.4,459.4000000000001,792])
saveas(gcf, './Outputs/Figure4_V5_part2_negcoop.svg')

% Make same figure for L-asp data
figure()
nrows = round(length(Files1)./2); % Keep same dimensions are for meAsp data
for jj = 1:length(Files2)
    %Load Data
    load([convertStringsToChars(Files2(jj)), 'plotData.mat'])

    Lplot = CDFPlotData.Lplot;

    %Generate sample curves for error bars
    rnd = CDFPlotData.MCMCSamples; %Dimensions: [mu, sigma, A (not used)]
    
    LplotPDF = log(min(CDFPlotData.Lplot)):0.01:log(max(CDFPlotData.Lplot));
    %Get error on pdft
     curveSamples = zeros(size(rnd,1), length(LplotPDF));
    for j = 1:size(rnd, 1)
        curveSamples(j, :) = PDFFunc(rnd(j, :), LplotPDF);
    end
    curvePosErr = prctile(curveSamples, 97.5);
    curveNegErr = prctile(curveSamples, 2.5);


    subplot(nrows, 2, jj)
    hold on
    y = PDFFunc(CDFPlotData.p_opt, LplotPDF);
%     plot(LplotPDF, y ./ max(y), 'k-', 'Linewidth', 2)
    shadeLineError(LplotPDF, curvePosErr, curveNegErr, [0.2, 0.2, 0.2])
    plot(LplotPDF, y, 'k-', 'Linewidth', 1)
%         ylabel(['normalized', newline, 'pdf'])
    xlim([log(xlimits(1)), log(xlimits(2))])
    ylabel('pdf')
%         if(jj == length(FilesTemp))
%             set(gca, 'Position', get(gca, 'Position') + [0, -0.02, 0, -0.1])
%         end
    meAsp_temp = meAsp_list2(jj); 
    Lasp_temp = Lasp_list2(jj);
    plot(LplotPDF, normpdf(LplotPDF, mu_Lasp(params, meAsp_temp, Lasp_temp), sig_Lasp(params, meAsp_temp, Lasp_temp)),...
        'r--', 'linewidth', 1)
end


set(gcf, 'Position', [750,-9.4,459.4000000000001,792])
saveas(gcf, './Outputs/Figure4_V5_part3_negcoop.svg')


% % %% Look at distribution of P0
% % P0_dist = @(x) lognpdf(x-1, p_opt(3), p_opt(4));
% % xrange = 1:0.001:3
% % 
% % figure();
% % plot(xrange, P0_dist(xrange))
% % set(gca, 'xscale', 'log')
% % 
% % P0_samps = 1 + lognrnd(p_opt(3), p_opt(4), [1000, 1]);
% % n_samps = log(exp(-0.8)+2)./log(P0_samps);
% % 
% % figure()
% % histogram(n_samps)
% 
% %% Calculate K0 and sigma0 from P0 distribution
% % Khalf functions
% Ki1 = 0.6238;
% Ki2 = 17.57;
% 
% f = @(L1, L2) (1 + L1./Ki1 + L2./Ki2);
% Khalf_Lasp = @(L1_0, L2_0, P0) (P0.*f(L1_0, L2_0) - L2_0./Ki2 - 1).*Ki1;
% Khalf_meAsp = @(L1_0, L2_0, P0) (P0.*f(L1_0, L2_0) - L1_0./Ki1 - 1).*Ki2;
% 
% 
% % Sample phenotypes
% nSamples = 10000;
% mu = -3.4285;
% sigma = 0.8037;
% rng(1)
% P0_samples = 1+lognrnd(mu, sigma, [nSamples, 1]);
% 
% % Calculate average and standard deviation of Khalf
% K0 = mean(Khalf_meAsp(0,0,P0_samples));
% sigma0 = std(Khalf_meAsp(0,0,P0_samples));








