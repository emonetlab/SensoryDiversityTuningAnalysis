%This function calculates response amplitudes for each measured response
%Response amplitude is defined as the average of post-stimulus FRET minus
%the average pre-stimulus FRET. 

%The current setup removes frames (badBefore and badAfter) relative to the
%first stimulus frame. 

%Updated 6/11/2021 JM to remove indices from the end of each measurement
%based on the parameter badEnd in parameters.m

function doseData = calcResponseAmp2(doseData, n, parms)
    for i = 1:length(doseData)
        preStimInd = doseData(i).s == doseData(i).s(1);
        postStimInd = doseData(i).s ~= doseData(i).s(1);
        
        %Modifiers to which indices to keep
        badBefore = parms.badBefore; %Bad indices before fist stim frame
        badAfter = parms.badAfter; %Bad indices after first stim frame
        badEnd = parms.badEnd; %Frames to remove from end of measurement
        
        AllInd = find(postStimInd);
        if ~(badBefore == 0 && badAfter == 0)
            badInd = [AllInd(1) - badBefore:AllInd(1) + badAfter];
            postStimInd(badInd) = 0;
        end
        
        if badEnd > 0
            badInd2 = [AllInd(end - badEnd + 1):AllInd(end)];
            postStimInd(badInd2) = 0;
        end
        
        
        %Allocate memory to store average pre and post-stimulus A
        doseData(i).preStimVals = zeros(size(doseData(i).a1, 1),n);
        doseData(i).postStimVals = zeros(size(doseData(i).a1, 1),n);
        
        %Calculate response amplitude as difference between mean pre and
        %post stimulus a
        %Rows are mean for each response. Columns are different stimulus
        %levels
%         doseData(i).preStimVals(:, 1) = mean(doseData(i).a1(:, preStimInd), 2);
%         doseData(i).postStimVals(:,1) = mean(doseData(i).a1(:, postStimInd), 2);
%         
%         if n > 1
%             doseData(i).preStimVals(:, 2) = mean(doseData(i).a2(:, preStimInd), 2);
%             doseData(i).postStimVals(:,2) = mean(doseData(i).a2(:, postStimInd), 2);
%         end
%         if n > 2
%             doseData(i).preStimVals(:, 3) = mean(doseData(i).a3(:, preStimInd), 2);
%             doseData(i).postStimVals(:,3) = mean(doseData(i).a3(:, postStimInd), 2);
%         end
%         if n > 3
%             doseData(i).preStimVals(:, 4) = mean(doseData(i).a4(:, preStimInd), 2);
%             doseData(i).postStimVals(:,4) = mean(doseData(i).a4(:, postStimInd), 2);
%         end
%         if n > 4
%             doseData(i).preStimVals(:, 5) = mean(doseData(i).a5(:, preStimInd), 2);
%             doseData(i).postStimVals(:,5) = mean(doseData(i).a5(:, postStimInd), 2);
%         end
       
        doseData(i).preStimVals(:, 1) = median(doseData(i).a1(:, preStimInd), 2);
        doseData(i).postStimVals(:,1) = median(doseData(i).a1(:, postStimInd), 2);
        
        if n > 1
            doseData(i).preStimVals(:, 2) = median(doseData(i).a2(:, preStimInd), 2);
            doseData(i).postStimVals(:,2) = median(doseData(i).a2(:, postStimInd), 2);
        end
        if n > 2
            doseData(i).preStimVals(:, 3) = median(doseData(i).a3(:, preStimInd), 2);
            doseData(i).postStimVals(:,3) = median(doseData(i).a3(:, postStimInd), 2);
        end
        if n > 3
            doseData(i).preStimVals(:, 4) = median(doseData(i).a4(:, preStimInd), 2);
            doseData(i).postStimVals(:,4) = median(doseData(i).a4(:, postStimInd), 2);
        end
        if n > 4
            doseData(i).preStimVals(:, 5) = median(doseData(i).a5(:, preStimInd), 2);
            doseData(i).postStimVals(:,5) = median(doseData(i).a5(:, postStimInd), 2);
        end

        %Raw response amplitude to be normalized later by dividing by
        %cell-scale a0
        doseData(i).rawR = doseData(i).postStimVals - doseData(i).preStimVals;
        
        %Response amplitude normalized by each measured response's
        %pre-stimulus a
        doseData(i).normR = repmat(median(doseData(i).postStimVals) ./ (median(doseData(i).preStimVals, "all")), [7,1]);
        doseData(i).a0 = median(doseData(i).preStimVals, "all");
%         median(doseData(i).normR)
%         doseData(i).normR = doseData(i).rawR ./ (median(doseData(i).preStimVals, "all")); 


    end
end