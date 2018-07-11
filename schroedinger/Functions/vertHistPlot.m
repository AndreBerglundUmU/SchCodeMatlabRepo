function vertHistPlot(storage,schemeShortNames,yAxisVector,leftTitle,axesInfo)
% vertHistPlot  - A function which will plot several histogram samples
%                 vertically, coupled with their titles.
% Syntax: vertHistPlot(storage,schemeShortNames,yAxisVector,leftTitle,axesInfo)
%
% Input:
% storage           - A Nxm matrix containing N samples of m random variables.
% schemeShortNames  - A cell vector of length m containing strings.
% yAxisVector       - A vector containing the y-axis labels.
% leftTitle         - A string containing the left side title.
% axesInfo          - A vector of length 4 containing the histogram sizes.
%                     Suggestion: [0.1305,0.774/m,0.11,0.815].
%
% Non-standard dependencies: None.
% See also: PSHist.m for example usage.

    % Plot a temporary figure which will help with the histogram axes
    tempFig = figure('visible','off');
    histogram(1,yAxisVector)
    ax = tempFig.CurrentAxes;
    pause(0.5)
    clear tempFig
    
    numSchemes = size(storage,2);
    figure
    % Should long labels be used, consider slanting them
    % xticklabels(temp(:,2))
    % xtickangle(45)
    
    % Set outer axis information and write out the labels
    axis([ 0 numSchemes ax.XLim])
    xticks((0:(numSchemes-1))+0.5)
    xticklabels(schemeShortNames)
        
    ylabel(leftTitle,'Interpreter','latex')
    set(gca,'FontSize',35)
    
    % Plot each histogram, removing the ticks and labels
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

