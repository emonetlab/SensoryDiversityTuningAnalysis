%This function combines the data in doseData in a form that is easier to
%plot. First, it combines each cells response amplitude data into a single
%3D table (allRData). Also combines all single-cell responses to each
%amplitude in a single nxm matrix, SCAvgData.allRData.

%It then calculates the average response amplitude at each level.

function [popAvgData, SCAvgData] = calcPopAvgDoseResponse2(m, doseData, parms, dataparms)
    %Get list of all concentration levels
    concList = unique(dataparms.concLevels);

    %Assemble all measured responses into one large 3D matrix
    popAvgData.allRData = NaN(m, length(concList), length(doseData));
    popAvgData.allRData_norm = NaN(m, length(concList), length(doseData));
    popAvgData.allPostStim = NaN(m, length(concList), length(doseData));
    for i = 1:length(doseData)
        popAvgData.allRData(:, find(ismember(concList,doseData(i).stimLevels)), i) = doseData(i).rawR; %(rep, doselevel, cell)
        popAvgData.allRData_norm(:, find(ismember(concList,doseData(i).stimLevels)), i) = doseData(i).normR; 
%         popAvgData.a0(i) = mean(doseData(i).preStimVals(:));
        popAvgData.a0(i) = median(doseData(i).preStimVals(:), 'all');
        popAvgData.allPostStim(:, find(ismember(concList,doseData(i).stimLevels)), i) = doseData(i).postStimVals;
    end
    
    %Assemble responses with all data in a single 2D matrix
    for i = 1:size(popAvgData.allRData, 2)
        temp = popAvgData.allRData(:, i, :);
        SCAvgData.allRData(:, i) = temp(:); %(rep, doselevel)
        
        temp = popAvgData.allRData_norm(:, i, :);
        SCAvgData.allRData_norm(:, i) = temp(:);
    end
    
    if ~parms.meanormedian
    	%Calculate mean response amplitude
        popAvgData.meanRraw = nanmean(popAvgData.allRData, [3, 1]);
        popAvgData.stdRraw = std(popAvgData.allRData, 0, [3, 1]);
        popAvgData.meanRnorm = nanmean(popAvgData.allRData_norm, [3, 1]);
        popAvgData.stdRnorm = std(popAvgData.allRData_norm, 0, [3, 1]);
        
        %Single-cell average dataset
        SCAvgData.ScMeanRraw = transpose(squeeze(nanmean(popAvgData.allRData, 1)));
        SCAvgData.ScStdRraw = transpose(squeeze(std(popAvgData.allRData,0, 1)));
        SCAvgData.ScMeanRnorm = transpose(squeeze(nanmean(popAvgData.allRData_norm, 1)));
        SCAvgData.ScStdRnorm = transpose(squeeze(std(popAvgData.allRData_norm,0, 1)));
    else
        %Calculate median response amplitude
        popAvgData.meanRraw = median(popAvgData.allRData, [3, 1]);
%         popAvgData.stdRraw = std(popAvgData.allRData, 0, [3, 1]);
        popAvgData.stdRraw = 1.4826*mad(popAvgData.allRData, 0, [3, 1]);
        popAvgData.meanRnorm = median(popAvgData.allRData_norm, [3, 1]);
%         popAvgData.stdRnorm = std(popAvgData.allRData_norm, 0, [3, 1]);
        popAvgData.stdRnorm = 1.4826*mad(popAvgData.allRData_norm, 0, [3, 1]);

        %Single-cell average dataset
        SCAvgData.ScMeanRraw = transpose(squeeze(median(popAvgData.allRData, 1)));
        SCAvgData.ScStdRraw = transpose(squeeze(std(popAvgData.allRData,0, 1)));
        SCAvgData.ScMeanRnorm = transpose(squeeze(median(popAvgData.allRData_norm, 1)));
        SCAvgData.ScStdRnorm = transpose(squeeze(std(popAvgData.allRData_norm,0, 1)));
    end
end