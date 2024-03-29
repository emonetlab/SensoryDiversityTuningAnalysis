%% Explanation:
%This script is to generate two plots for Figure S1
% Panel E. Simulated dose-response curves with response and K1/2
% distributions to show how we can infer one from the other
% Panel F. Example response to show how the response amplitude is
% calculated

%% Load data
clear
close all
dataparms = datasets_ForPublication();
% data = dataparms(1); %0 background meAsp responses
data = dataparms(2); %1 uM Ser background meAsp responses
% data = dataparms(3); %1 uM Ser background meAsp responses
addpath('./Functions')
addpath('./Data')
addpath('./Outputs')
parms = parameters(); %Struct containing system parameters
Files = dataparms.Files;
concLevels = dataparms.concLevels;
EfretData = combineFretFiles(Files, concLevels); %Combine files into one dataset



%% Figure 1a: Diagram of stimulus protocol
data1a = EfretData(1);
tall = [data1a.t_min_max(:); data1a.t(:)];
sall = [data1a.s_min_max(:); data1a.s(:)];

%Rename concentrations to correct levels
sall(sall == -1) = 100;
sall(sall == 0) = 100;
sall(sall == 2) = 102;
sall(sall == 10) = 105;
sall(sall == 20) = 110;
sall(sall == 40) = 120;
sall(sall == 80) = 140;
sall(sall == 1000) = 200;

%Sort vectors to make good stairs plot
[tsorted, torder] = sort(tall);
ssorted = sall(torder);
for i = 1:length(ssorted)-1
    if ssorted(i+1) == 100
        ssorted(i) = 100;
    end
end

% figure();
% sp = subplot(1, 3, 2)
% %Plot to scale with figure 1b
% stairs(tsorted/60, ssorted + 20, 'Linewidth', 1.5, 'Color', [0.4, 0.4, 0.4])
% % set(gca, 'yscale', 'log')
% ylabel('[L]')
% % xlabel('Time (min)')
% % xlim([-1, 32])
% xlim([3, 28.5])
% ylim([100, 210])
% sp.Position(1) = 0.16;
%     sp.Position(3) = 0.63;
% set(gca,'YTick',[], 'XTick', [])
% set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 16)
% % set(gcf, 'Position', [489 343 450 350])
% set(gcf, 'Position', [155.4,312.2,1228,50])
% % JMAxes


%% Figure 1b: Plot multiple cell's full time-series right on top of each other
selectedCells = [43, 2, 74, 11, 5];
ncells = length(selectedCells);
nPanels = 3*ncells;
extralength = .5;

%Get plot-friendly data for stimulus time-series
plotData = getPlotFriendlyData(EfretData);
switchTimes = find(plotData.sChanges);
avgTSeries = mean(plotData.tSeries, 1)./60;
switchTimes_sat1 = find(plotData.s_sat1Changes);
avgTSat1 = mean(plotData.t_sat1Series, 1)./60;
switchTimes_sat2 = find(plotData.s_sat2Changes);
avgTSat2 = mean(plotData.t_sat2Series, 1)./60;

f = figure();
for ii = 1:ncells
    %Plot left saturating stimulus
    sp =subplot(ncells, 3, (ii-1)*3 + 1);
    LplotDataT = EfretData(selectedCells(ii)).t_min_max(1, :)/60;
    LplotDataA = EfretData(selectedCells(ii)).smooth_a_min_max(1,:);
    plot(LplotDataT, LplotDataA, '.', 'Color', '#e02459')
    xlim([min(LplotDataT), max(LplotDataT)]) 
    ylim([-0.3, 1])
    %Color stimulus regions
    for jj = 1:2:length(switchTimes_sat1)
        if jj < length(switchTimes_sat1)
            x = avgTSat1(transpose(switchTimes_sat1(jj:jj+1) - [0;2])); %Had to subtract 2 because of nans
        else
            x = [avgTSat1(switchTimes_sat1(jj)), max(avgTSat1)];
        end
        
        y = get(gca,'YLim');
        nx = [x(1) x(2) x(2) x(1)];
        ny = [y(1) y(1) y(2) y(2)];
        hold on
        patch(nx,ny, 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'FaceColor', [.5, .5, .5]);
    end
%     sp.Position = sp.Position + [extralength + 0.02 0 -0.17 0];
    sp.Position(3) = 0.02;
    if ii ~= ncells
        set(gca,'XTick',[])
    end
    if ii == round(ncells/2)
        ylabel('a')
    end
    set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 12)

    %     Plot right saturating stimulus
    sp = subplot(ncells, 3, (ii-1)*3 + 3);
    LplotDataT = EfretData(selectedCells(ii)).t_min_max(2, :)/60;
    LplotDataA = EfretData(selectedCells(ii)).smooth_a_min_max(2,:);
    plot(LplotDataT, LplotDataA, '.', 'Color', '#e02459')
    xlim([min(LplotDataT), max(LplotDataT)])  
    ylim([-0.3, 1])
%     Color stimulus regions
    for jj = 1:2:length(switchTimes_sat2)
        if jj < length(switchTimes_sat2)
            x = avgTSat2(transpose(switchTimes_sat2(jj:jj+1) - [0;2])); %Had to subtract 2 because of nans
        else
            x = [avgTSat2(switchTimes_sat2(jj)), max(avgTSat2)];
        end
        
        y = get(gca,'YLim');
        nx = [x(1) x(2) x(2) x(1)];
        ny = [y(1) y(1) y(2) y(2)];
        hold on
        patch(nx,ny, 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'FaceColor', [.5, .5, .5]);
    end
%     sp.Position = sp.Position + [-0.15 0 -0.17 0];
%     sp.Position = sp.Position + [extralength 0 -0.17 0];
    sp.Position(3) = 0.02;
    sp.Position(1) = 0.8;
    set(gca,'YTick',[])
    if ii ~= ncells
        set(gca,'XTick',[])
    end
    set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 12)

    
    %Plot center time-series
    sp = subplot(ncells, 3, (ii-1)*3 + 2);
    CplotDataT = EfretData(selectedCells(ii)).t/60;
    CplotDataA = EfretData(selectedCells(ii)).smooth_a;   
    plot(CplotDataT', CplotDataA', '-', 'Linewidth', 1.5, 'Color', '#e02459')
    xlim([min(CplotDataT(:)), max(CplotDataT(:))]) 
    ylim([-0.3, 1])
    %Color stimulus regions
    for jj = 1:2:length(switchTimes)
        if jj < length(switchTimes)
            x = avgTSeries(transpose(switchTimes(jj:jj+1) - [0;2])); %Had to subtract 2 because of nans
        else
            x = [avgTSeries(switchTimes(jj)), max(avgTSeries)];
        end
        
        y = get(gca,'YLim');
        nx = [x(1) x(2)+2/60 x(2)+2/60 x(1)];
        ny = [y(1) y(1) y(2) y(2)];
        hold on
        patch(nx,ny, 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', [.6, .6, .6]);
    end
%     sp.Position = sp.Position + [-0.15 0 extralength 0];
%     sp.Position = sp.Position + [0 0 extralength 0];
    sp.Position(1) = 0.16;
    sp.Position(3) = 0.63;
    set(gca,'YTick',[])
    if ii ~= ncells
        set(gca,'XTick',[])
    else
        xlabel('Time (min)')
    end
    set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 12)
    
    
end
% set(gcf, 'Position', [489 343 800 350])
% f.Position = f.Position + [0 -100 0 100];

set(gcf, 'Position', [155.4,312.2,1228,408.8])

%% Figure 1b+ plot first cell zoomed in somewhere in the time-series
selectedCell = selectedCells(4);

figure();
CplotDataT = EfretData(selectedCell).t/60;
CplotDataA = EfretData(selectedCell).smooth_a;
CplotDataA_real = EfretData(selectedCell).a;

% sp = subplot(2, 1, 1);
% stairs(tsorted/60, ssorted, 'Linewidth', 1.5)
% % set(gca, 'yscale', 'log')
% ylabel('[L]')
% xlim([14, 17.5])
% ylim([100, 150])
% set(gca,'YTick',[], 'XTick',[])
% sp.Position = sp.Position + [0, -.1, 0, -.2];
% set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 16)
% 
% subplot(2, 1, 2)
hold on
plot(CplotDataT', CplotDataA_real', 'o', 'Color', [0.2, 0.2, 0.2])
% xlim([min(CplotDataT(:)), max(CplotDataT(:))]) 
ylim([-0.3, .8])
%Color stimulus regions
for jj = 1:2:length(switchTimes)
    if jj < length(switchTimes)
        x = avgTSeries(transpose(switchTimes(jj:jj+1) - [0;2])); %Had to subtract 2 because of nans
    else
        x = [avgTSeries(switchTimes(jj)), max(avgTSeries)];
    end

    y = get(gca,'YLim');
    nx = [x(1) x(2)+2/60 x(2)+2/60 x(1)];
    ny = [y(1) y(1) y(2) y(2)];
    hold on
    patch(nx,ny, 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', [.8, .8, .8]);
end
plot(CplotDataT', CplotDataA_real', 'o', 'Color', [0.2, 0.2, 0.2])
plot(CplotDataT', CplotDataA', '-', 'Linewidth', 1.5, 'Color', '#e02459')
xlabel("Time (min)")
ylabel("a")
xlim([15.5, 16.8])
set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 12)
set(gcf, 'Position', [583.4,454.6,465.6,200])

%% Figure 1d: Single-cell representative average response
nLevels = 5;
nResponses = 7;
selectedCell = 22;
respNum = 5;

cellPlotData = EfretData(selectedCell);

selectedResponses = zeros(nResponses, size(cellPlotData.a, 2));
for ii = 1:nResponses
    selectedResponses(ii, :) = cellPlotData.a(nLevels.*(ii-1)+respNum, :);
end
figure(); hold on
plot(0:0.5:9.5, (selectedResponses), 'o', 'Color', [0.8, 0.8, 0.8])
plot(0:0.5:9.5, (median(selectedResponses, 1)), 'Linewidth', 3, 'Color', '#e02459')
yline(median(selectedResponses(:, 1:11),'all'), 'linewidth', 2)
yline(median(selectedResponses(:, 12:end), 'all'), 'linewidth', 2)
% plot(0:0.5:9.5, median(selectedResponses, 1), 'o', 'Color', '#e02459')
% plot(0:0.5:9.5, (selectedResponses), 'o', 'Color', '#e02459')
xlabel("t (s)")
ylabel("a")
set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 12)
set(gcf, 'Position', [554,338,192,168], 'Renderer', 'painters')
saveas(gcf, './Outputs/FigureS1F.svg')
saveas(gcf, './Outputs/FigureS1F.png')


%% Figure 1c: Explanation of K1/2 CDF from single-cell responses

%Constants
Ki = 18; %uM
Ka = 2900; %uM
alpha = 2; %kT
m0 = 0.5; 
muK = 7;
sigmaK = 2;


%N and mstar Distribution 
Ksamples = lognrnd(muK, sigmaK, 30, 1);
mstar = lognrnd(.01, .001, 30, 1);

%Function for kinase activity
activity = @(L, K) 1./(1 + L./K);

%For each value of n, calculate dose-response curve
L = 10.^[0:0.001:6];


figure()
for ii = 1:length(Ksamples)
    subplot(2, 2, 3)
    hold on
    A = activity(L, Ksamples(ii));
    Atest = activity(10^3, Ksamples(ii));
    if (max(A)-Atest)./(max(A)-min(A)) >= 0.5
        col = '#e02459';
    else 
        col = '#c2c2c2';
    end
    plot(L, activity(L, Ksamples(ii)), 'Color', col , 'Linewidth', 1);
end
xline(10^3, '-.', 'Linewidth', 2)
yline(0.5, '-.', 'Linewidth', 2)
set(gca, 'xscale', 'log')
set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 16)
set(gca, 'XTick', [])
xlabel("[L_i]")
ylabel("R")

sp = subplot(2, 2, 1)
L1 = 10.^[0:0.001:3];
hold on
area(L1, lognpdf(L1, muK + 3, sigmaK), 'FaceColor','#e02459' )
plot(L, lognpdf(L, muK + 3, sigmaK),'k', 'Linewidth', 2)
set(gca, 'xscale', 'log')
xline(10^3, '-.', 'Linewidth', 2)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel("P(K_{1/2})")
set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 14)
sp.Position = sp.Position + [0 -.08 0 -.2];

sp = subplot(2, 2, 4)
Pl = 0:0.01:1;
Pl1 = .5:0.01:1;
hold on
area(Pl1, normpdf(Pl1, 0.5, 0.15), 'FaceColor','#e02459' )
plot(Pl, normpdf(Pl, 0.5, 0.15), 'k', 'Linewidth', 2)
xline(0.5, '-.', 'Linewidth', 2)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
sp.Position = sp.Position + [-.08 0 -0.2 0];
ylabel("P(R(L_i))")
set(gca,'TickLabelInterpreter', 'latex', 'FontName', 'Arial', 'FontSize', 14)

camroll(-90)


set(gcf, 'Position', [248.2,162.6,690.8000000000001,530.4])


