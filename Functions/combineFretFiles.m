%This script does 1 thing:
    %It takes in a list of FRET data files, and combines them into a single
    %file so data can be analyzed all at once.

function EfretDataFinal = combineFretFiles(Files, concLevels)
    %Load first EfretData file
    load(Files(1), "reorgData");
    EfretDataFinal = reorgData.resp_data;
    %Assign each cell a list of concentrations that were measured in it
    %Assign each cell a number based on which file it originated from
    for j = 1:length(EfretDataFinal)
        EfretDataFinal(j).concLevels = concLevels(1, :);
        EfretDataFinal(j).DataSet = 1;
    end

    for i = 2:length(Files)
        load(Files(i), "reorgData");
        EfretDataTemp = reorgData.resp_data;
        
        %Add conc levels to new cell's data
        for j = 1:length(EfretDataTemp)
            EfretDataTemp(j).concLevels = concLevels(i, :);
            EfretDataTemp(j).DataSet = i;
        end
        
        %Update all fields in EfretData with new cell's data
        for j = 1:length(EfretDataTemp)
            EfretDataFinal(length(EfretDataFinal) + 1) = EfretDataTemp(j);
        end      
    end
end