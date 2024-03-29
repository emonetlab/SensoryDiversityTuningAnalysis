%% This makes figure 5 to describe simulation results
%Will begin with a diagram of the simulation environment, and possible a
%diagram of the n distribution sampled to generate cells. Then, will show
%the performance of a population with diverse n compared to the median
%phenotype.

close all
clear
clc

% SimData = {'./Outputs/SimulationData/NobackPerformanceData.mat',...
%     './Outputs/SimulationData/UncompPerformanceData.mat',...
%     './Outputs/SimulationData/CompPerformanceData.mat'};

% SimData = {'./Outputs/SimulationData/NobackPerformanceData.mat'};
% SimData = {'./Outputs/SimulationData/PerformanceData2Gradient.mat'};


%% Plot fraction of population at each peak for diverse and non-diverse
%First, do 2 uncompetitive gradients.
SimData = {'./Outputs/SimulationData/PerformanceDataUnComp220726.mat'};
% colr = {"#10FFA2",...
%     "#10E5FF",...
%     "#FF10E5"};

colr = {"#4E4E4E",...
    "#10E5FF",...
    "#FF10E5"};

meAspStrengthChoice = 100;

figure()

for ii = 1:length(SimData)
    load(SimData{ii})

    subplot(4, 2, 1)
    %Choose only one meAsp source strength
    sourcePlot = SourceStrengthNoDiv_sort(:, 1) == meAspStrengthChoice;
    
    %Plotting parameters
    lw = 3;
    
    %Make plot

    hold on
    plot(SourceStrengthNoDiv_sort(sourcePlot, 2), smooth(NoDiversity_sort(sourcePlot, 1))./10000, '--', 'linewidth', lw, 'Color', "#058ED9")
    plot(SourceStrengthNoDiv_sort(sourcePlot, 2), smooth(NoDiversity_sort(sourcePlot, 2))./10000, '--', 'linewidth', lw, 'Color', "#FF934F")

    plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(Diversity_sort(sourcePlot, 1))./10000, '-', 'linewidth', lw, 'Color', "#058ED9")
    plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(Diversity_sort(sourcePlot, 2))./10000, '-', 'linewidth', lw, 'Color', "#FF934F")

    legend(["meAsp Peak", "serine peak"], 'Location','best')
    set(gca,'xscale', 'log', 'yscale', 'log', ...
        'XTick', [10^-3, 10^-2, 10^-1, 10^0, 10^1], ...
        'YTick', [0.02, 0.05, 0.1, 0.2, 0.4])
    xlabel("Serine Source strength \muM", 'fontsize', 14)
    ylabel("Fraction at source", 'fontsize', 14)
    ylim([10^-2, 0.5])


   %Make plot of average n_Tar value that makes it to each peak
    avgNTarAtPeak = zeros(size(SourceStrengthNoDiv_sort, 1), 2); %[Concentration, peak no.]
    for jj = 1:size(SourceStrengthNoDiv_sort, 1)
        cellsInPeak1_ind = cellsInCenter_Div_sort(:, jj, 1);
        cellsInPeak2_ind = cellsInCenter_Div_sort(:, jj, 2);
        nValuesAtPeak1 = cellNValues_Div_sort(logical(cellsInPeak1_ind), jj, 1);
        nValuesAtPeak2 = cellNValues_Div_sort(logical(cellsInPeak2_ind), jj, 1);
        avgNTarAtPeak(jj, 1) = median(nValuesAtPeak1);
        avgNTarAtPeak(jj, 2) = median(nValuesAtPeak2);
    end
    
    subplot(4, 2, 3)
    plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(avgNTarAtPeak(sourcePlot, 1)), 'Color', "#058ED9", 'linewidth', lw)
    hold on
    plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(avgNTarAtPeak(sourcePlot, 2)), 'Color', "#FF934F", 'linewidth', lw)
    set(gca,'xscale', 'log', 'yscale', 'log', ...
        'XTick', [10^-3, 10^-2, 10^-1, 10^0, 10^1])
    ylabel("\langle n_{Tar} \rangle")

    %Make plot of average n_Tsr value that makes it to each peak
    avgNTarAtPeak = zeros(size(SourceStrengthNoDiv_sort, 1), 2); %[Concentration, peak no.]
    for jj = 1:size(SourceStrengthNoDiv_sort, 1)
        cellsInPeak1_ind = cellsInCenter_Div_sort(:, jj, 1);
        cellsInPeak2_ind = cellsInCenter_Div_sort(:, jj, 2);
        nValuesAtPeak1 = cellNValues_Div_sort(logical(cellsInPeak1_ind), jj, 2);
        nValuesAtPeak2 = cellNValues_Div_sort(logical(cellsInPeak2_ind), jj, 2);
        avgNTarAtPeak(jj, 1) = median(nValuesAtPeak1);
        avgNTarAtPeak(jj, 2) = median(nValuesAtPeak2);
    end
    
    subplot(4, 2, 5)
    plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(avgNTarAtPeak(sourcePlot, 1)), 'Color', "#058ED9", 'linewidth', lw)
    hold on
    plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(avgNTarAtPeak(sourcePlot, 2)), 'Color', "#FF934F", 'linewidth', lw)
    set(gca,'xscale', 'log', 'yscale', 'log', ...
        'XTick', [10^-3, 10^-2, 10^-1, 10^0, 10^1])
    ylabel("\langle n_{Tsr} \rangle")


    % Ratio of diverse to non-diverse population
    subplot(4, 2, 7)
    Ratio1 = smooth(Diversity_sort(sourcePlot, 1))./smooth(NoDiversity_sort(sourcePlot, 1));
    Ratio2 = smooth(Diversity_sort(sourcePlot, 2))./smooth(NoDiversity_sort(sourcePlot, 2));
    hold on
    plot(SourceStrengthNoDiv_sort(sourcePlot, 2), Ratio1, 'linewidth', lw, 'Color', "#058ED9")
    plot(SourceStrengthNoDiv_sort(sourcePlot, 2), Ratio2, 'linewidth', lw, 'Color', "#FF934F")
    ylim([0, 2.2])
    ylabel("Relative performance")
    set(gca,'xscale', 'log', ...
        'XTick', [10^-3, 10^-2, 10^-1, 10^0, 10^1])


end


%% Next do 2 competitive gradients
SimData = {'./Outputs/SimulationData/PerformanceData220726.mat'};
% colr = {"#10FFA2",...
%     "#10E5FF",...
%     "#FF10E5"};

colr = {"#4E4E4E",...
    "#10E5FF",...
    "#FF10E5"};

meAspStrengthChoice = 100;


ii = 1;
load(SimData{ii})

subplot(4, 2, 2)
%Choose only one meAsp source strength
sourcePlot = SourceStrengthNoDiv_sort(:, 1) == meAspStrengthChoice;

%Plotting parameters
lw = 3;

%Make plot of performance as a function of L-asp strength
hold on
plot(SourceStrengthNoDiv_sort(sourcePlot, 2), smooth(NoDiversity_sort(sourcePlot, 1))./10000, '--', 'linewidth', lw, 'Color', "#058ED9")
plot(SourceStrengthNoDiv_sort(sourcePlot, 2), smooth(NoDiversity_sort(sourcePlot, 2))./10000, '--', 'linewidth', lw, 'Color', "#e02459")

plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(Diversity_sort(sourcePlot, 1))./10000, '-', 'linewidth', lw, 'Color', "#058ED9")
plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(Diversity_sort(sourcePlot, 2))./10000, '-', 'linewidth', lw, 'Color', "#e02459")

legend(["meAsp Peak", "L-Asp peak"], 'Location','best')
set(gca,'xscale', 'log', 'yscale', 'log', ...
    'XTick', [10^-3, 10^-2, 10^-1, 10^0, 10^1], ...
    'YTick', [0.02, 0.05, 0.1, 0.2, 0.4])
xlabel("L-Asp Source strength \muM", 'fontsize', 14)
ylabel("Fraction at source", 'fontsize', 14)
ylim([10^-2, 0.5])

%Make plot of average n_Tar value that makes it to each peak
avgNTarAtPeak = zeros(size(SourceStrengthNoDiv_sort, 1), 2); %[Concentration, peak no.]
for jj = 1:size(SourceStrengthNoDiv_sort, 1)
    cellsInPeak1_ind = cellsInCenter_Div_sort(:, jj, 1);
    cellsInPeak2_ind = cellsInCenter_Div_sort(:, jj, 2);
    nValuesAtPeak1 = cellNValues_Div_sort(logical(cellsInPeak1_ind), jj, 1);
    nValuesAtPeak2 = cellNValues_Div_sort(logical(cellsInPeak2_ind), jj, 1);
    avgNTarAtPeak(jj, 1) = median(nValuesAtPeak1);
    avgNTarAtPeak(jj, 2) = median(nValuesAtPeak2);
end

subplot(4, 2, 4)
plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(avgNTarAtPeak(sourcePlot, 1)), 'Color', "#058ED9", 'linewidth', lw)
hold on
plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(avgNTarAtPeak(sourcePlot, 2)), 'Color', "#e02459", 'linewidth', lw)
set(gca,'xscale', 'log', 'yscale', 'log', ...
    'XTick', [10^-3, 10^-2, 10^-1, 10^0, 10^1])
ylabel("\langle n_{Tar} \rangle")

%Make plot of average n_Tsr value that makes it to each peak
avgNTarAtPeak = zeros(size(SourceStrengthNoDiv_sort, 1), 2); %[Concentration, peak no.]
for jj = 1:size(SourceStrengthNoDiv_sort, 1)
    cellsInPeak1_ind = cellsInCenter_Div_sort(:, jj, 1);
    cellsInPeak2_ind = cellsInCenter_Div_sort(:, jj, 2);
    nValuesAtPeak1 = cellNValues_Div_sort(logical(cellsInPeak1_ind), jj, 2);
    nValuesAtPeak2 = cellNValues_Div_sort(logical(cellsInPeak2_ind), jj, 2);
    avgNTarAtPeak(jj, 1) = median(nValuesAtPeak1);
    avgNTarAtPeak(jj, 2) = median(nValuesAtPeak2);
end

subplot(4, 2, 6)
plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(avgNTarAtPeak(sourcePlot, 1)), 'Color', "#058ED9", 'linewidth', lw)
hold on
plot(SourceStrengthDiv_sort(sourcePlot, 2), smooth(avgNTarAtPeak(sourcePlot, 2)), 'Color', "#e02459", 'linewidth', lw)
set(gca,'xscale', 'log', 'yscale', 'log', ...
    'XTick', [10^-3, 10^-2, 10^-1, 10^0, 10^1])
ylabel("\langle n_{Tsr} \rangle")

 % Ratio of diverse to non-diverse population
subplot(4, 2, 8)
Ratio1 = smooth(Diversity_sort(sourcePlot, 1))./smooth(NoDiversity_sort(sourcePlot, 1));
Ratio2 = smooth(Diversity_sort(sourcePlot, 2))./smooth(NoDiversity_sort(sourcePlot, 2));
hold on
plot(SourceStrengthNoDiv_sort(sourcePlot, 2), Ratio1, 'linewidth', lw, 'Color', "#058ED9")
plot(SourceStrengthNoDiv_sort(sourcePlot, 2), Ratio2, 'linewidth', lw, 'Color', "#e02459")
set(gca,'xscale', 'log', ...
    'XTick', [10^-3, 10^-2, 10^-1, 10^0, 10^1])
ylim([0, 2.2])

set(gcf, 'Position', [510,45,785,835])
saveas(gcf, './Outputs/Figure6a.svg')

%% Plot environment
figure()
resolution = 1000;
x = linspace(-500, 500, resolution);
y = linspace(-500, 500, resolution);
[X, Y] = meshgrid(x, y);
Z = 100 .* exp(-10^(-2.75) .* ((X-100).^2 + (Y-100).^2));
Z2 = 100 .* exp(-10^(-2.75) .* ((X+100).^2 + (Y+100).^2));
contour(X, Y, Z-100, 'linewidth', 2)
hold on
contour(X, Y, Z2+50, 'linewidth', 2)
xline(0)
yline(0)
xlim([-200, 200])
ylim([-200, 200])
xlabel("X (\mum)")
ylabel("Y (\mum)")


%Load example cell trajectories
AgentData = './Outputs/SimulationData/ExampleTrajectories220725.mat';
load(AgentData)
x_pos = squeeze(posOut(:, 1, 1:100));
y_pos = squeeze(posOut(:, 2, 1:100));

cellno = [7, 2, 9];
for i = 1:length(cellno)
    plot(x_pos(cellno(i), :), y_pos(cellno(i), :), 'linewidth', .8)
end

set(gcf, 'Position', [29,342,660,551])
saveas(gcf, './Outputs/Figure6b.svg')




