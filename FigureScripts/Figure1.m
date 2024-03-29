%% Explanation:
%This script is to generate plots for the first figure of my diversity
%tuning paper. The key plot is meant to display the average response of
%several individual cells to 5 different levels. This will take the form of
%an n_cell x 5 grid of plots, where each plot contains the 7 measurements
%for that cell at that concentration shown in gray, and the average
%response shown in bold. color to be chosen later.
clear
close all

%% Plot parameters
colors = {'#FF934F', '#058ED9', 'black'};


%=======================
%% 0 Background Version
%=======================

%% Load data into a friendly form
dataparms = datasets_ForPublication();
data = dataparms(1); %0 background meAsp responses
% data = dataparms(3); %1 uM Ser background meAsp responses
addpath('./Functions')
addpath('./Data')
addpath('./Outputs')
parms = parameters(); %Struct containing system parameters
Files = data.Files;
concLevels = data.concLevels;
EfretData = combineFretFiles(Files, concLevels); %Combine files into one dataset


% For each cell, make a new matrix: (nReps : Time : nChannels)
% This will contain all the responses per stimulus level for that cell
nChannels = 5;
nReps = size(EfretData(1).a, 1) ./ nChannels;
nPoints = size(EfretData(1).a, 2);
ChannelInd = reshape(repmat(1:nChannels, nReps, 1)', [], 1);

for ii = 1:length(EfretData)
    EfretData(ii).a_bychannel = zeros(nReps, nPoints, nChannels);
    for channel = 1:nChannels
        EfretData(ii).a_bychannel(:, :, channel) = EfretData(ii).a(ChannelInd == channel,:);
%         EfretData(ii).a_bychannel(:, :, channel) = EfretData(ii).smooth_a(ChannelInd == channel,:);
    end  
end


%% Plot grid of responses for select cells
figure()
% selectedCells = [43, 5, 1, 11, 5];
selectedCells = [5, 19, 1, 43, 11];
gridRows = length(selectedCells) + 1;
%Make a grid that is nCells X nStimulusLevels

%Plot stimulus in the top row
s = EfretData(1).s;
s(s==2) = 0.2;
s(s==10) = 0.5;
s(s==20) = 1;
s(s==40) = 2;
s(s==80) = 4;

for stimlevel = 1:nChannels
    %Plot stimulus L1
    subplot(gridRows, nChannels, stimlevel)
    stairs(s(stimlevel, :), 'color', colors{2}, 'linewidth', 1.5)
    ylim([0, 6])
    set(gca, 'YTick', [0, 3, 6], 'XTick', [], 'XColor', 'none')
    box off

    
    %Plot stimulus L2
    subplot(gridRows, nChannels, stimlevel)
    hold on
    yline(0, 'linewidth', 1.5, 'color', colors{1})
%     ylim([0,2])
%     set(gca, 'YTick', [0, 2], 'XTick', [], 'XColor', 'none')
    
    %Handle axis labels
    if stimlevel == 1
        ylabel("[L]")
    else
        set(gca, 'YTick', [])
    end
end



%Plot measurements
for row = 1:length(selectedCells)
    celNo = selectedCells(row);
    for stimlevel = 1:nChannels
        %Plot FRET measurements
        subplot(gridRows, nChannels, (row)*nChannels + stimlevel)
        a_temp = EfretData(celNo).a_bychannel;
        plot(a_temp(:, :, stimlevel)', 'linewidth', 1, 'color', [0.8, 0.8, 0.8])
        hold on
        plot(mean(a_temp(:, :, stimlevel), 1), 'color', colors{3}, 'linewidth', 1.5)
        ylim([-0.3, 1])

        %Handle axis labels
        if row < length(selectedCells)
            set(gca, 'XTick', [], 'XColor', 'none')
        end
        
        if stimlevel > 1
            set(gca, 'YTick', [])
        end

        if stimlevel == 3 && row == length(selectedCells) 
            xlabel("Time (sec)")
        elseif stimlevel == 1 && row == round(length(selectedCells)/2)
            ylabel('Activity')
        end
        
        box off
    end
end
set(gcf, 'Position', [533,262,620,648],'Renderer', 'painters')

saveas(gcf, './Outputs/Figure1_1.svg')



%==============================
%% Non-zero Background Version
%==============================
dataparms = datasets_ForPublication();
% data = dataparms(1); %0 background meAsp responses
data = dataparms(3); %1 uM Ser background meAsp responses
addpath('./Functions')
addpath('./Data')
addpath('./Outputs')
parms = parameters(); %Struct containing system parameters
Files = data.Files;
concLevels = data.concLevels;
EfretData = combineFretFiles(Files, concLevels); %Combine files into one dataset


% For each cell, make a new matrix: (nReps : Time : nChannels)
% This will contain all the responses per stimulus level for that cell
nChannels = 5;
nReps = size(EfretData(1).a, 1) ./ nChannels;
nPoints = size(EfretData(1).a, 2);
ChannelInd = reshape(repmat(1:nChannels, nReps, 1)', [], 1);

for ii = 1:length(EfretData)
    EfretData(ii).a_bychannel = zeros(nReps, nPoints, nChannels);
    for channel = 1:nChannels
        EfretData(ii).a_bychannel(:, :, channel) = EfretData(ii).a(ChannelInd == channel,:);
%         EfretData(ii).a_bychannel(:, :, channel) = EfretData(ii).smooth_a(ChannelInd == channel,:);
    end  
end


%% Plot grid of responses for select cells
figure()
selectedCells = [1, 4, 9, 10, 12];
gridRows = length(selectedCells) + 1;
%Make a grid that is nCells X nStimulusLevels

%Plot stimulus in the top row
s = EfretData(1).s;
s(s==2) = 0.2;
s(s==10) = 0.5;
s(s==20) = 1;
s(s==40) = 2;
s(s==80) = 4;

for stimlevel = 1:nChannels
    %Plot stimulus L1
    subplot(gridRows, nChannels, stimlevel)
    stairs(s(stimlevel, :), 'color', colors{2}, 'linewidth', 1.5)
    ylim([0, 6])
    set(gca, 'YTick', [0, 3, 6], 'XTick', [], 'XColor', 'none')
    box off

    
    %Plot stimulus L2
    subplot(gridRows, nChannels, stimlevel)
    hold on
    yline(1, 'linewidth', 1.5, 'color', colors{1})
%     ylim([0,2])
    
    %Handle axis labels
    if stimlevel == 1
        ylabel("[L]")
    else
        set(gca, 'YTick', [])
    end
end



%Plot measurements
for row = 1:length(selectedCells)
    celNo = selectedCells(row);
    for stimlevel = 1:nChannels
        %Plot FRET measurements
        subplot(gridRows, nChannels, (row)*nChannels + stimlevel)
        a_temp = EfretData(celNo).a_bychannel;
        plot(a_temp(:, :, stimlevel)', 'linewidth', 1, 'color', [0.8, 0.8, 0.8])
        hold on
        plot(mean(a_temp(:, :, stimlevel), 1), 'color', colors{3}, 'linewidth', 1.5)
        ylim([-0.3, 1])

        %Handle axis labels
        if row < length(selectedCells)
            set(gca, 'XTick', [], 'XColor', 'none')
        end
        
        if stimlevel > 1
            set(gca, 'YTick', [])
        end

        if stimlevel == 3 && row == length(selectedCells) 
            xlabel("Time (sec)")
        elseif stimlevel == 1 && row == round(length(selectedCells)/2)
            ylabel('Activity')
        end
        
        box off
    end
end
set(gcf, 'Position', [533,262,620,648], 'Renderer', 'painters')

saveas(gcf, './Outputs/Figure1_2.svg')
