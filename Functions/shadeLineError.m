function shadeLineError(Lplot, curvePosErr, curveNegErr, lineColor)
    %This function take in an x axis (Lplot) and the upper and lower bounds
    %of the error for a line, then plots a shaded area between the two
    %curves.
%     plot(Lplot, curveNegErr, 'linewidth', .1, 'Color', [0,0,0])
%     plot(Lplot, curvePosErr, 'linewidth', .1, 'Color', [0,0,0])
    x2 = [Lplot, fliplr(Lplot)];
    inBetween = [curveNegErr, fliplr(curvePosErr)];
    fill(x2, inBetween, 'g', 'faceColor', lineColor, 'facealpha',.5, 'LineStyle','none')
end