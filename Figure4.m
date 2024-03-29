close all; clear;
addpath('./Functions')

figure()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ====Fitting Khalf CV Values=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Piece 1: plot the CV curves.
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
    "./Outputs/1LAsp_1MeAspBack_MeAspDR/",...
    "./Outputs/1LAsp_10MeAspBack_MeAspDR/",...
    "./Outputs/100GluBackMeAspDR/",...
    "./Outputs/100MeAsp_1uMMaltose_BackMeAspDR/"];

% Gave 0 background small value for plot, but fit with 0
bckg = [0, 10, 100, 0.3, 0.1, 1, 0.01, 10, 1, 100, 1, 10, 0, 0];
bckg_plot = bckg; bckg_plot(bckg==0) = 10^-3;
LAspBckg = [0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 2, 2, -1, -1];

% Extract CV from each file, and store in order in vectors
% [CV, means, stds, err] = extract_cv_from_files(CVFiles);
[CV, means, stds, CVerr, CVlowerbound, CVupperbound, meanerr,...
    meanlowbound, meanupperbound] = extract_cv_from_files(CVFiles);


subplot(4, 2, [1, 3])
hold on
%Parameters
% KiMeAsp = 18; %uM
% KiLAsp = 0.5; %uM
KiMeAsp = 17.556; %uM
KiLAsp = 0.624; %uM


% nmu = 2.82;
% nsigma = 1.15;
% nTar = lognrnd(nmu,nsigma, 1000, 1);

%Regime I assumption
%Assuming the cluster is in regime I, K1/2 = (ln2)/n .* Ki
%When the background increases, Kieff = Ki(1 + L0/Ki)
% Khalf = @(n, L0, L02, Ki1, Ki2) log(2)./n .* Ki1 .* (1 + L0./Ki1 + L02./Ki2);
% FCKhalf = @(n, L0, L02, Ki1, Ki2) (Khalf(n, L0, L02, Ki1, Ki2)) ./ L0;

%Calculate FCKhalf curve as a function of background
colr = ["#4E4E4E"; "#798086"; "#50B2C0"];
mrkrSize = 7;
% 
L0Range = 10.^[-3:0.1:3];
% KhalfRange0 = Khalf(nTar, L0Range, 0, KiMeAsp, KiLAsp) + L0Range;
% KhalfRange1 = Khalf(nTar, L0Range, 1, KiMeAsp, KiLAsp)+ L0Range;
% KhalfRange10 = Khalf(nTar, L0Range, 10, KiMeAsp, KiLAsp)+ L0Range;
% 
% CV0 = std(KhalfRange0) ./ mean(KhalfRange0);
% CV1 = std(KhalfRange1) ./ mean(KhalfRange1);
% CV10 = std(KhalfRange10) ./ mean(KhalfRange10);

K0 = means(1);
sigma0 = stds(1); 
disp('meAsp experimental parameters:')
K0
sigma0
% sigma0 = means(1) * 1.5;

% K0 = 1.44;
% sigma0 = 2.418;
% CVfunc = @(L0, Lasp)  sigma0 ./ (K0 + L0./(1 + L0./KiMeAsp + Lasp./KiLAsp));
% CVfunc = @(L0, Lasp)  sigma0 ./ (K0 + L0./(1 + L0./7.5 + Lasp./0.75));

%p = [sigma0, K0, Ki_measp)
% CVfitfunc = @(p, L0) p(1) ./ (p(2) + L0./(1 + L0./p(3)));
% CVfitfunc_competitor = @(p, L0, Lasp) p(1) ./ (p(2) + L0./(1 + L0./p(3) + Lasp./p(4)));

% Fix values of Ki to literature values
% CVfitfunc = @(p, L0) p(1) ./ (p(2) + L0./(1 + L0./KiMeAsp));
% CVfitfunc_competitor = @(p, L0, Lasp) p(1) ./ (p(2) + L0./(1 + L0./KiMeAsp + Lasp./KiLAsp));

% Fix the value of K0, since we have high confidence in this
% CVfitfunc = @(p, L0) p(1) ./ (K0 + L0./(1 + L0./p(3))); 
% CVfitfunc_competitor = @(p, L0, Lasp) p(1) ./ (K0 + L0./(1 + L0./p(3) + Lasp./p(4)));

% Fix the value of K0 and Ki
CVfitfunc = @(p, L0) p(1) ./ (K0 + L0./(1 + L0./KiMeAsp)); 
CVfitfunc_competitor = @(p, L0, Lasp) p(1) ./ (K0 + L0./(1 + L0./KiMeAsp + Lasp./KiLAsp));

backgrounds = bckg(LAspBckg==0); backgrounds2 = bckg(LAspBckg==2); backgrounds10 = bckg(LAspBckg==10);
y_fit = CV(LAspBckg==0); y_fit2 = CV(LAspBckg==2); y_fit10 = CV(LAspBckg==10);
err_fit = CVerr(LAspBckg==0); err_fit2 = CVerr(LAspBckg==2); err_fit10 = CVerr(LAspBckg==10);

p0 = [sigma0, K0, KiMeAsp, KiLAsp]; %Initial parameters
% Log posterior including error at all 3 background L-asp concentrations
inv_log_post = @(p) sum((y_fit - CVfitfunc(p, backgrounds)).^2./err_fit.^2./2) +... 
    sum((y_fit2 - CVfitfunc_competitor(p, backgrounds2, 2)).^2./(2*err_fit2.^2)) + ...
    sum((y_fit10 - CVfitfunc_competitor(p, backgrounds10, 10)).^2./(2*err_fit10.^2));
% inv_log_post = @(p) sum((y_fit - CVfitfunc(p, backgrounds)).^2./err_fit.^2./2);
p_opt = fminunc(inv_log_post, p0);
% print('Optimal CV Curve Parameters')
disp('MeAsp Response Parameters: [sigma0, K0, KimeAsp, KiLasp]')
p0
p_opt % [sigma0, K0, Ki]
meAsp_CV_params = p_opt;

%===Plot with free K0 and Sigma0
CV0 = CVfitfunc(p_opt, L0Range);
CV10 = CVfitfunc_competitor(p_opt, L0Range, 10);
CV2 = CVfitfunc_competitor(p_opt, L0Range, 2);
% CV0 = CVfitfunc(p0, L0Range);
% CV10 = CVfitfunc_competitor(p0, L0Range, 10);
% CV2 = CVfitfunc_competitor(p0, L0Range, 2);

hold on
C(1) = plot(L0Range, CV0, 'Linewidth', 2, 'Color', colr(1,:));
C(2) = plot(L0Range, CV10, 'Linewidth', 2, 'Color', colr(3,:));
C(3) = plot(L0Range, CV2, 'Linewidth', 2, 'Color', colr(2,:));
errorbar(bckg_plot(LAspBckg==0), CV(LAspBckg==0),  CV(LAspBckg==0)-CVlowerbound(LAspBckg==0),...
    CVupperbound(LAspBckg==0)-CV(LAspBckg==0), 'o', 'MarkerEdgeColor', colr(1,:),...
    'Color', colr(1,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
errorbar(bckg_plot(LAspBckg==10), CV(LAspBckg==10), CV(LAspBckg==10)-CVlowerbound(LAspBckg==10),...
    CVupperbound(LAspBckg==10)-CV(LAspBckg==10), 'o', 'MarkerEdgeColor', colr(3,:),...
    'Color', colr(3,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
errorbar(bckg_plot(LAspBckg==2), CV(LAspBckg==2), CV(LAspBckg==2)-CVlowerbound(LAspBckg==2),...
    CVupperbound(LAspBckg==2)-CV(LAspBckg==2), 'o', 'MarkerEdgeColor', colr(2,:),...
    'Color', colr(2,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
ylabel("Predicted K_{1/2} CV")
xlabel("Background [meAsp] \muM")
xlim([0.001, 1000])
% ylim([0,.75])
xticks([10^-2, 10^-1, 10^-0, 10^1, 10^2])
% set(gca, 'xscale', 'log', ...
%     'Position', get(gca, 'Position') + [0, +0.03, 0, -0.15])
set(gca, 'xscale', 'log')
hold off

%% Plot Mean and Std of K1/2
subplot(4, 2, 2)
hold on
plot(bckg(LAspBckg==0), means(LAspBckg==0), 'o', 'Color', colr(1,:))
ylabel("Mean K_{1/2}")
xlabel("[meAsp] \muM")
xlim([-10, 110])
ylim([-10, 110])
% xlim([0.001, 1000])
% ylim([0,.75])
% xticks([10^-2, 10^-1, 10^-0, 10^1, 10^2])
% set(gca, 'xscale', 'log', 'yscale', 'log')

subplot(4, 2, 4)
hold on
plot(bckg(LAspBckg==0), stds(LAspBckg==0), 'o', 'Color', colr(1,:))
ylabel("Std Deviation in K_{1/2}")
xlabel("[meAsp] \muM")
xlim([-10, 110])
% ylim([-1, 7])
% xlim([0.001, 1000])
% ylim([0,.75])
% xticks([10^-2, 10^-1, 10^-0, 10^1, 10^2])
% set(gca, 'xscale', 'log', 'yscale', 'log')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ====LAsp Responses=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVFiles = ["./Outputs/0BackLAspDose/", ...
    "./Outputs/1LAsp_LAspDR/", ...
    "./Outputs/0_01LAsp_LAspDR/", ...
    "./Outputs/0_5LAsp_LAspDR/",...
    "./Outputs/0_1LAsp_LAspDR/", ...
    "./Outputs/0_05LAsp_LAspDR/",...
    "./Outputs/10LAsp_LAspDR/",...
    "./Outputs/1LAsp_100meAsp_LAspDR/",...
    "./Outputs/10LAsp_100meAsp_LAspDR_V2/",...
    "./Outputs/100MeAspBackLAspDose/",...
    "./Outputs/100meAsp_0_1LAsp_LAspDR/"];

bckg = [0, 1, 0.01, 0.5, 0.1, 0.05, 10, 1, 10, 0, 0.1];
bckg_plot = bckg; bckg_plot(bckg==0) = 10^-3;
meAspBckg = [0, 0, 0, 0, 0, 0, 0, 100, 100, 100, 100];

% [CV, means, stds, err] = extract_cv_from_files(CVFiles);
[CV, means, stds, CVerr, CVlowerbound, CVupperbound, meanerr,...
    meanlowbound, meanupperbound] = extract_cv_from_files(CVFiles);

subplot(4, 2, [5, 7])
hold on
%Parameters
KiMeAsp = meAsp_CV_params(3); %uM
KiLAsp = meAsp_CV_params(4); %uM

colr = ["#4E4E4E"; "#848FA2"; "#848FA2"];
L0Range = 10.^[-3:0.1:3];

%p = [sigma0, K0, Ki_measp)
K0 = means(1);
sigma0 = CV(1).*K0; 

% % All parameters free
% CVfitfunc = @(p, L0) p(1) ./ (p(2) + L0./(1 + L0./p(3)));
% CVfitfunc_competitor = @(p, L0, Measp) p(1) ./ (p(2) + L0./(1 + L0./p(3) + Measp./p(4)));

% Fix K0
% CVfitfunc = @(p, L0) p(1) ./ (K0 + L0./(1 + L0./p(3)));
% CVfitfunc_competitor = @(p, L0, Measp) p(1) ./ (K0 + L0./(1 + L0./p(3) + Measp./p(4)));

% Fix K0 and Ki
CVfitfunc = @(p, L0) p(1) ./ (K0 + L0./(1 + L0./KiLAsp));
CVfitfunc_competitor = @(p, L0, Measp) p(1) ./ (K0 + L0./(1 + L0./KiLAsp + Measp./KiMeAsp));

% % Fix Kimeasp
% CVfitfunc_fixedK = @(p, L0) p(1) ./ (p(2) + L0./(1 + L0./KiLAsp));
% CVfitfunc_competitor_fixedK = @(p, L0, Measp) p(1) ./ (p(2) + L0./(1 + L0./p(3) + Measp./KiMeAsp));

% Fix KimeAsp and K0
CVfitfunc_fixedK = @(p, L0) p(1) ./ (K0 + L0./(1 + L0./KiLAsp));
CVfitfunc_competitor_fixedK = @(p, L0, Measp) p(1) ./ (K0 + L0./(1 + L0./p(3) + Measp./KiMeAsp));

backgrounds = bckg(meAspBckg==0);
y_fit = CV(meAspBckg==0);
err_fit = CVerr(meAspBckg==0);
backgrounds2 = bckg(meAspBckg==100);
y_fit2 = CV(meAspBckg==100);
err_fit2 = CVerr(meAspBckg==100);
p0 = [sigma0, K0, KiLAsp, KiMeAsp];
% With unfixed K values
inv_log_post = @(p) sum((y_fit - CVfitfunc(p, backgrounds)).^2./(2*err_fit.^2)) +...
    1*sum((y_fit2 - CVfitfunc_competitor(p, backgrounds2, 100)).^2./(2*err_fit2.^2));
p_opt = fminunc(inv_log_post, p0);
disp('LAsp Response Parameters: [sigma0, K0, KiLasp, KiMeAsp]')
p0
p_opt % [sigma0, K0, Ki]

% With fixed K values
inv_log_post = @(p) sum((y_fit - CVfitfunc_fixedK(p, backgrounds)).^2./(2*err_fit.^2)) +...
    0*sum((y_fit2 - CVfitfunc_competitor_fixedK(p, backgrounds2, 100)).^2./(2*err_fit2.^2));
p_opt_fixedK = fminunc(inv_log_post, p0);

%===Plot with free K0 and Sigma0
CV0 = CVfitfunc(p_opt, L0Range);
CV100 = CVfitfunc_competitor(p_opt, L0Range, 100);
CV100_fixedK = CVfitfunc_competitor_fixedK(p_opt_fixedK, L0Range, 100);
C(1) = plot(L0Range, CV0, 'Linewidth', 2, 'Color', colr(1,:));
C(2) = plot(L0Range, CV100, 'Linewidth', 2, 'Color', colr(2,:));
C(3) = plot(L0Range, CV100_fixedK, '--', 'Linewidth', 2, 'Color', colr(2,:));
errorbar(bckg_plot(meAspBckg==0), CV(meAspBckg==0), CV(meAspBckg==0)-CVlowerbound(meAspBckg==0),...
    CVupperbound(meAspBckg==0)-CV(meAspBckg==0), 'o', 'MarkerFaceColor', colr(1,:), 'MarkerEdgeColor', [0.9 0.9 0.9],...
    'Color', colr(1,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
errorbar(bckg_plot(meAspBckg==100), CV(meAspBckg==100), CV(meAspBckg==100)-CVlowerbound(meAspBckg==100),...
    CVupperbound(meAspBckg==100)-CV(meAspBckg==100), 'o', 'MarkerFaceColor', colr(2,:), 'MarkerEdgeColor',  [0.9 0.9 0.9],...
    'Color', colr(2,:), 'Linewidth', 0.9, 'MarkerSize', mrkrSize)
ylabel("Predicted K_{1/2} CV")
xlabel("Background [L-Asp] \muM")
xlim([0.001, 1000])
% ylim([0,.75])
xticks([10^-2, 10^-1, 10^-0, 10^1, 10^2])
% set(gca, 'xscale', 'log', ...
%     'Position', get(gca, 'Position') + [0, +0.03, 0, -0.15])
ylim([0, 1.7])
set(gca, 'xscale', 'log')
hold off

%% Plot Mean and Std of K1/2
subplot(4, 2, 6)
hold on
plot(bckg(meAspBckg==0), means(meAspBckg==0), 'o', 'Color', colr(1,:))
ylabel("Mean K_{1/2}")
xlabel("[L-Asp] \muM")
xlim([-1, 11])
ylim([-1, 12])
% xlim([0.001, 1000])
% ylim([0,.75])
% xticks([10^-2, 10^-1, 10^-0, 10^1, 10^2])
% set(gca, 'xscale', 'log', 'yscale', 'log')

subplot(4, 2, 8)
hold on
plot(bckg(meAspBckg==0), stds(meAspBckg==0), 'o', 'Color', colr(1,:))
ylabel("Std Deviation in K_{1/2}")
xlabel("[L-Asp] \muM")
xlim([-1, 11])
ylim([-.1, 0.8])
% xlim([0.001, 1000])
% ylim([0,.75])
% xticks([10^-2, 10^-1, 10^-0, 10^1, 10^2])
% set(gca, 'xscale', 'log', 'yscale', 'log')

%% Figure aesthetics
% set(gcf, 'Position', [151.4,97,896.6,665])
set(gcf, 'Position', [151.4,2,896.6,750])

saveas(gcf, './Outputs/Figure4_V5_part1.svg')
saveas(gcf, './Outputs/Figure4_V5_part1.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting meAsp khalf distribution with Lasp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Files are organized by row:
    %Row 1: top left
    %Row 2: bottom left
    %Row 3: top right
    %Row 4: bottom right
% Files1 = {["./Outputs/MeAspDR_OneConcentrationSets/", ...
%     "./Outputs/10MeAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/", ...
%     "./Outputs/10LAsp_100meAsp_MeAspDR/"]};
Files1 = {["./Outputs/MeAspDR_OneConcentrationSets/", ...
    "./Outputs/10MeAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/"]};
% Files1 = {["./Outputs/0BackLAspDose/", ...
%     "./Outputs/10MeAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/", "./Outputs/10MeAsp10LAspBackMeAspDR/"]};



% Color Palette To Color Each plot. Rows organized the same way
colrs = {["#4E4E4E", "#058ED9", "#e02459", "#058ED9", "#058ED9"]};

lnstyles = {["-", "-", "-", "--", "--", "--"]};

showMarkers = {[1, 1, 1, 0, 1]};

backgrColors = {[1 1 1]};

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
        mu = CDFPlotData.p_opt(1);
        sigma = CDFPlotData.p_opt(2);

        m_all = exp(rnd(:, 1) + rnd(:,2).^2./2);
        v_all = exp(2*rnd(:, 1) + rnd(:, 2).^2).*(exp(rnd(:, 2).^2) - 1);

    end
end
%%
% Plot predicted khalf distributions using fitted meAsp model parameters
CVfitfunc_competitor = @(p, L0, Lasp) p(1) ./ (p(2) + L0./(1 + L0./p(3) + Lasp./p(4)));
Meanfunc_competitor = @(p, L0, Lasp) p(2)*(1+L0./p(3) + Lasp./p(4)) + L0;
params = meAsp_CV_params; %[sigma0, K0, Ki_measp, Ki_lasp]
% params = [1.14, 1.2, 64, 1.1];
% params(3) = 3;

% Get lognormal parameters from mean and CV
m = @(params, L0, Lasp) Meanfunc_competitor(params, L0, Lasp);
v = @(params, L0, Lasp) (m(params, L0, Lasp).*CVfitfunc_competitor(params, L0, Lasp)).^2;
mu = @(params, L0, Lasp) log(m(params, L0, Lasp).^2./(sqrt(v(params, L0, Lasp)+m(params, L0, Lasp).^2)));
sig = @(params, L0, Lasp) sqrt(log(v(params, L0, Lasp)./m(params, L0, Lasp).^2 + 1));
% 
% Lrange = 10.^[-3:0.01:3];
% subplot(4, 2, plotInds(1, 2))
% hold on
% L0 = 0; Lasp = 0;
% plot(Lrange, logncdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp)),  'k--', 'linewidth', 2)
% L0 = 10; Lasp = 0;
% plot(Lrange, logncdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp)),  'k--', 'linewidth', 2)
% L0 = 10; Lasp = 10;
% plot(Lrange, logncdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp)), 'k--', 'linewidth', 2)
% L0 = 100; Lasp = 10;
% plot(Lrange, logncdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp)), 'k--', 'linewidth', 2)
% set(gca, 'xscale', 'log', 'XTick', [10^-2, 10^-1, 10^0, 10^1, 10^2])
% xlim([10.^-3, 10.^3])
% 
subplot(4, 2, plotInds(1, 1))
% hold on
% Lrange = -4:0.01:4;
% L0 = 0; Lasp = 0;
% plot(Lrange, normpdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp))./max(normpdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp))),  'k--', 'linewidth', 2)
% 
% Lrange = 2:0.01:6;
% L0 = 10; Lasp = 0;
% plot(Lrange, normpdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp))./max(normpdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp))),  'k--', 'linewidth', 2)
% 
% Lrange = 2:0.01:6;
% L0 = 10; Lasp = 10;
% plot(Lrange, normpdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp))./max(normpdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp))),  'k--', 'linewidth', 2)
% 
% 
% Lrange = 3:0.01:7;
% L0 = 100; Lasp = 10;
% plot(Lrange, normpdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp))./max(normpdf(Lrange, mu(params, L0, Lasp), sig(params, L0, Lasp))), 'k--', 'linewidth', 2)
set(gca, 'Position', get(gca, 'Position') + [0, -0.02, 0, -0.1])
% 

set(gcf, 'Position', [1011.4,137.8,422.4000000000001,478.4])
saveas(gcf, './Outputs/Figure4_V5_part2.svg')
saveas(gcf, './Outputs/Figure4_V5_part2.png')
