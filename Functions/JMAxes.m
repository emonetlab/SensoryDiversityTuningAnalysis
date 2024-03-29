% JMAxes.m
% makes current figure pretty: based on PrettyFig.m
% 1. makes all line widths of plots 2

function [] = JMAxes(varargin)
    warning off
    % defaults
    lw = 2; % line width of graphical elements
    plw = 2; % plot line width 
    fs = 16; % font size

    % get handle to all plots in current figure
    axesHandles = findall(gcf,'type','axes');
    
    for i = 1:length(axesHandles)
        % set line width and font size
%         set(axesHandles(i),'FontSize',fs,'LineWidth',lw, 'TickLabelInterpreter', 'latex', 'FontName', 'Arial')
        set(axesHandles(i),'FontSize',fs, 'TickLabelInterpreter', 'latex', 'FontName', 'Arial')

        % turn the minor ticks on
        set(axesHandles(i),'XMinorTick','on','YMinorTick','on')	
        
        % find all plots and set those line widths appropriately
        ph=get(axesHandles(i),'Children');
        for j = 1:length(ph)
            try
                set(ph(j),'LineWidth',plw)
            catch
                % probably an image or something.
                % so reverse tick direction
                set(gca,'TickDir','out')
                box on
            end
        end
        
        % Add all axis tiks to log scales
        if  length(get(axesHandles(i),'XTick')) <= 3 && strcmp(get(axesHandles(i),'XScale'),'log')
            c=get(axesHandles(i),'Children');
            minlog = Inf;maxlog = -Inf;
            for k = 1:length(c)
                minlog = min([ min(nonzeros(get(c(k),'XData'))) minlog]);
                maxlog = max([ max(nonzeros(get(c(k),'XData'))) maxlog]);
            end
            a = ceil(log10(minlog));
            z = floor(log10(maxlog));
            if length(a:z) > 2
                set(axesHandles(i),'XTick',10.^(a:z));
            end
        end
    end
    
   
end