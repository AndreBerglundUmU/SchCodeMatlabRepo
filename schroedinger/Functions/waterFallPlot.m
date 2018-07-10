function waterFallPlot(varargin)
    sideFun = @(u) absSq(u);
    if nargin == 5
        XInt = varargin{2};
        K = varargin{3};
        xVec = linspace(XInt(1),XInt(2),K);
    elseif nargin > 5
        xVec = varargin{6};
        if nargin == 7
            sideFun = @(u) abs(u);
        end
    end
    run = varargin{1};
    T = varargin{4};
    N = varargin{5};
    tVec = linspace(0,T,N);
    temp = waterfall(xVec,tVec,sideFun(run));
    set(temp,'LineWidth',1.5)
    ylabel('t','FontSize',20);
    xlabel('x','FontSize',20);
    if nargin > 5
        zlabel('$|u_n|$','Rotation',0,'HorizontalAlignment','right','Interpreter','latex','FontSize',20);
    else
        zlabel('$|u_n|^2$','Rotation',0,'HorizontalAlignment','right','Interpreter','latex','FontSize',20);
    end
    % str=sprintf('Space-time evolution. $h=$ %.2e, $M=$ %d',dt,Nx);
    % title(str,'Interpreter','latex','FontSize',14)
    axis tight
end