%This script takes EFRET data, and using the number of replicate
%measurements of each stimulus level n, and the number of stimulus levels
%m, makes a struct, doseData, that contains each cell's responses for each stimulus
%level in their own table.

function doseData = reorganizeData(EfretData, n, m, OutputDest)
    nRows = size(EfretData(1).a, 1);
    nReps = nRows ./ n;
    sortMat = repmat(transpose(1:n), nReps, 1);
    
    %For each cell, reformat data
    for i = 1:size(EfretData, 2)
        %Get the stimulus protocol for that cell
        doseData(i).s = EfretData(i).s(2, :);
        %Get the stimulus levels used for that cell
        doseData(i).stimLevels = EfretData(i).concLevels;
        %Set up the first cell's response measurements
        doseData(i).a1 = EfretData(i).a(sortMat == 1, :);
        doseData(i).t1 = EfretData(i).t(sortMat == 1, :);
        doseData(i).a1 = doseData(i).a1(1:m, :);
        doseData(i).t1 = doseData(i).t1(1:m, :);
        aTemp(i, :, 1) = mean(doseData(i).a1, 1);
        
        if n > 1
            doseData(i).a2 = EfretData(i).a(sortMat == 2, :); 
            doseData(i).t2 = EfretData(i).t(sortMat == 2, :);
            doseData(i).a2 = doseData(i).a2(1:m, :);
            doseData(i).t2 = doseData(i).t2(1:m, :);
            aTemp(i, :, 2) = mean(doseData(i).a2, 1);
        end
        
        if n > 2
            doseData(i).a3 = EfretData(i).a(sortMat == 3, :); 
            doseData(i).t3 = EfretData(i).t(sortMat == 3, :);
            doseData(i).a3 = doseData(i).a3(1:m, :);
            doseData(i).t3 = doseData(i).t3(1:m, :);
            aTemp(i, :, 3) = mean(doseData(i).a3, 1);
        end
        
        if n > 3
            doseData(i).a4 = EfretData(i).a(sortMat == 4, :); 
            doseData(i).t4 = EfretData(i).t(sortMat == 4, :);
            doseData(i).a4 = doseData(i).a4(1:m, :);
            doseData(i).t4 = doseData(i).t4(1:m, :);
            aTemp(i, :, 4) = mean(doseData(i).a4, 1);
        end
        
        if n > 4
            doseData(i).a5 = EfretData(i).a(sortMat == 5, :); 
            doseData(i).t5 = EfretData(i).t(sortMat == 5, :);
            doseData(i).a5 = doseData(i).a5(1:m, :);
            doseData(i).t5 = doseData(i).t5(1:m, :);
            aTemp(i, :, 5) = mean(doseData(i).a5, 1);
        end
        
        if n>5
            doseData(i).a6 = EfretData(i).a(sortMat == 6, :); 
            doseData(i).t6 = EfretData(i).t(sortMat == 6, :);
            doseData(i).a6 = doseData(i).a6(1:m, :);
            doseData(i).t6 = doseData(i).t6(1:m, :);
            aTemp(i, :, 6) = mean(doseData(i).a6, 1);
        end
            
    end
    
    %% Plot average response to each stimlus level
    for i = 1:length(EfretData)
        DataSets(i) = EfretData(i).DataSet;
    end
    
    for j = 1:length(unique(DataSets))
        nPlots = size(aTemp, 3);
        figure()
        for i = 1:nPlots
    %         subplot(nPlots, 1, i)
            hold on
            plot(mean(aTemp(DataSets==j,:,i), 1), 'Linewidth', 3)
    %         title(['Channel ', num2str(i)])
            legend(["1", "2", "3", "4", "5", "6"], 'Location', 'southwest')
        end
        xlabel("Frames")
        ylabel("a")
        title("Population-average Responses")
        savefig(gcf, [OutputDest, 'AverageResponse', num2str(j),'.fig'])
        saveas(gcf, [OutputDest,  'AverageResponse', num2str(j), '.png'])
    end

end
