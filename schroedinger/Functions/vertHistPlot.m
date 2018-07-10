function vertHistPlot(storage,schemeShortNames,yAxisVector,leftTitle,axesInfo)

    % Plot a temporary figure which will help with the histogram axes
    tempFig = figure('visible','off');
    histogram(1,yAxisVector)
    ax = tempFig.CurrentAxes;
    pause(0.5)
    clear tempFig
    
    numSchemes = size(storage,2);
    figure
    %xticklabels(temp(:,2))
    % xtickangle(45)
    axis([ 0 numSchemes ax.XLim])
    xticks((0:(numSchemes-1))+0.5)
    xticklabels(schemeShortNames)
        
    %title(figTitle,'Interpreter','latex')
    ylabel(leftTitle,'Interpreter','latex')
    set(gca,'FontSize',35)

    for i = 1:numSchemes
        axes('Position',[axesInfo(1)+(i-1)*axesInfo(2) axesInfo(3) axesInfo(2) axesInfo(4)])
        histogram(storage(:,i),yAxisVector)
        set(gca,'view',[90 -90])
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        set(gca,'YTickLabel',[])
        set(gca,'XTickLabel',[])
        pause(0.1)
    end
end

