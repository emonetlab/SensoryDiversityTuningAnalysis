function parms = parameters()
    %% Show plots or not. Plots are still saved in output directory
    parms.showPlot = 'on'; %Changing to 'off' will suppress figures. 
    
    %% Analysis Parameters
    parms.m = 7; %Number of responses to include in analysis
    parms.meanormedian = 1; %If 0, use mean for population average. If 1, use median 
    parms.nSamples = 1000; %Number of samples for MCMC sampling
    parms.thining = 10;
    parms.nBootstraps = 1000; %Number of bootstrap samples for fitting CDF
    
    %% Removed frames
    %Determines how many frames before and after the first stimulus frame
    %are removed when calculating response amplitude
    parms.badBefore = 0;
    parms.badAfter = 1;
%     parms.badEnd = 0;
%     parms.badEnd = 2;
    parms.badEnd = 0;
%     parms.badEnd = 5;
    
    %% Plot parameters
    %Color blind friendly pallette #E1DAAE, #FF934F, #CC2D35, #058ED9, #848FA2, #2D3142
    parms.lineColor = "#CC2D35";
    parms.pointColor = "#058ED9";
    
end