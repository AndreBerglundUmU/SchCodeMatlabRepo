function waterFallPlot(run,XInt,TInt)
% waterFallPlot  -  A function which will plot a waterfall plot of the 
%                   absolute value squared.
% Syntax: waterFallPlot(run,XInt,TInt)
%
% Input:
% run   - A Nxm matrix containing N samples of m random variables.
% XInt  - A vector containing the end points of the spacial interval.
% TInt  - A vector containing the end points of the time interval.
%
% Non-standard dependencies: absSq.m.
% See also: FDWaterfall.m
%           PSWaterfall.m
%           for example usage.

    N = size(run,1);
    M = size(run,2);
    
    t = linspace(TInt(1),TInt(2),N);
    x = linspace(XInt(1),XInt(2),M);
    
    temp = waterfall(x,t,absSq(run));
    set(temp,'LineWidth',1.5)
    ylabel('t','FontSize',20);
    xlabel('x','FontSize',20);
    zlabel('$|u_n|^2$','Rotation',0,'HorizontalAlignment','right','Interpreter','latex','FontSize',20);
    % str=sprintf('Space-time evolution. $h=$ %.2e, $M=$ %d',dt,Nx);
    % title(str,'Interpreter','latex','FontSize',14)
    axis tight
end