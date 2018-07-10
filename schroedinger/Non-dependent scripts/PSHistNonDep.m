%% Set initial info and function
initSeed = 1;
batchSize = 1000;
% Time and area
% L = 30; XInt = [-L,L];
XInt = [0,2*pi];
T = 1/8; TInt = [0,T];

% Numerical precision, number of points
N = 2^10; % Time
M = 2^10; % Space

h = (TInt(2)-TInt(1))/N;
dx = (XInt(2)-XInt(1))/M;

alpha = 1;
q = 1;
c = 1;
u0Fun = @(x) sqrt(2*alpha/q)*exp(0.5*1i*c*x).*sech(sqrt(alpha)*x);
% u0Fun = @(x) 1.4*exp(-4*x.^2);

% Nonlinearity information
sigma = 1;
G = @(u) absSq(u).^sigma.*u;

% See trapezoidalIntegral.m
L2norm = @(currU) trapezoidalIntegral(dx,absSq([currU currU(1)]));

k = 2*pi/(XInt(2)-XInt(1))*[0:M/2-1, 0, -M/2+1:-1];
kSq = k.^2;

t = linspace(TInt(1),TInt(2),N+1);
x = XInt(1) + dx*(0:M-1);
u0FunVal = u0Fun(x);
u0L2Val = L2norm(u0FunVal);

% See makePSSchroedSchemes.m for schemes
numAvailableSchemes = 9;
schemesUsed = false(numAvailableSchemes,1);
schemesUsed(1) = true; % FEul
schemesUsed(2) = true; % BEul
schemesUsed(3) = true; % MEul
schemesUsed(4) = true; % CN
schemesUsed(5) = true; % EExp
schemesUsed(6) = true; % SExp
schemesUsed(7) = true; % LTSpl
% schemesUsed(8) = true; % FSpl
schemesUsed(9) = true; % SSpl

schemes.shortNames = cell(9,1);
schemes.shortNames{1} = 'FEul';
schemes.shortNames{2} = 'BEul';
schemes.shortNames{3} = 'MEul';
schemes.shortNames{4} = 'CN';
schemes.shortNames{5} = 'EExp';
schemes.shortNames{6} = 'SExp';
schemes.shortNames{7} = 'LTSpl';
schemes.shortNames{8} = 'FSpl';
schemes.shortNames{9} = 'SSpl';

numUsedSchemes = sum(schemesUsed);
schemeIndexMat = [(1:numUsedSchemes)' , find(schemesUsed)];

%% Query storage
maxL2DriftBatch = zeros(batchSize,numUsedSchemes);

%% Perform calculations
rng(initSeed,'twister')
parfor m = 1:batchSize
    % Load the broadcast variables to internal
    internalSchemeIndexMat = schemeIndexMat;
    
    W = randn(N,2)*sqrt(h/2);
    
    maxL2Drift = zeros(1,numUsedSchemes);
    for j = 1:numUsedSchemes
        currU = fft(u0FunVal);
        currScheme = internalSchemeIndexMat(j,2);
        for i = 1:N
            dW = W(i,:);
            %% Scheme and query calculations
            % See makePSSchroedSchemes.m for schemes
            switch currScheme
                case 1 % Forwards/explicit Euler scheme
                    realSpaceCurrU = ifft(currU);
                    a = 1i*sum(dW)/2*kSq;
                    b = (1 - a)./(1 + a).*currU;
                    c = 1i*h./(1 + a);
                    currU = b + c.*fft(G(realSpaceCurrU));
                    
                case 2 % Backwards/implicit Euler scheme
                    % Implicit solving for next currU
                    crit = true;
                    nextU = currU;
                    realSpaceCurrU = ifft(currU);
                    realSpaceNextU = realSpaceCurrU;
                    
                    a = 1i*sum(dW)/2*kSq;
                    b = (1 - a)./(1 + a).*currU;
                    c = 1i*h./(1 + a);
                    
                    noRounds = 1;
                    while crit && noRounds<120
                        tempU = nextU;
                        nextU = b + c.*fft(G(realSpaceNextU));
                        realSpaceNextU=ifft(nextU);
                        crit = norm((tempU-nextU)./M,2) > eps;
                        noRounds = noRounds+1;
                    end
                    currU = nextU;
                case 3 % Midpoint scheme
                    % Implicit solving for next currU
                    crit = true;
                    nextU = currU;
                    realSpaceCurrU = ifft(currU);
                    realSpaceNextU = realSpaceCurrU;
                    
                    a = 1i*sum(dW)/2*kSq;
                    b = (1 - a)./(1 + a).*currU;
                    c = 1i*h./(1 + a);
                    
                    noRounds = 1;
                    while crit && noRounds<120
                        tempU = nextU;
                        nextU = b + c.*fft(G((realSpaceNextU+realSpaceCurrU)/2));
                        realSpaceNextU=ifft(nextU);
                        crit = norm((tempU-nextU)./M,2) > eps;
                        noRounds = noRounds+1;
                    end
                    currU = nextU;
                case 4 % Crank-Nicolson scheme
                    % Implicit solving for next currU
                    crit = true;
                    nextU = currU;
                    realSpaceCurrU = ifft(currU);
                    realSpaceNextU = realSpaceCurrU;
                    
                    a = 1i*sum(dW)/2*kSq;
                    b = (1 - a)./(1 + a).*currU;
                    c = 1i*h./(1 + a);
                    
                    noRounds = 1;
                    while crit && noRounds<120
                        tempU = nextU;
                        nextU = b + c.*fft(CNnonLin(realSpaceNextU,realSpaceCurrU,sigma));
                        realSpaceNextU=ifft(nextU);
                        crit = norm((tempU-nextU)./M,2) > eps;
                        noRounds = noRounds+1;
                    end
                    currU = nextU;
                case 5 % Explitic exponential scheme
                    currU = exp(-sum(dW)*1i*kSq).*(currU + 1i*h*fft(G(ifft(currU))));
                case 6 % Symmetric exponential scheme
                    % Implicit solving for NStar
                    a = exp(-dW(1)*1i*kSq).*currU;
                    crit = true;
                    NStar = currU;
                    noRounds = 1;
                    while crit && noRounds<120
                        oldNStar = NStar;
                        tempNStar = a+h/2*NStar;
                        NStar = fft(1i*G(ifft(tempNStar)));
                        crit = norm((oldNStar-NStar)./M,2) > eps;
                        noRounds = noRounds+1;
                    end
                    % Last step
                    currU = exp(-sum(dW)*1i*kSq).*currU + h*exp(-dW(2)*1i*kSq).*NStar;
                case 7 % Lie-Trotter splitting scheme
                    % Full linear step
                    tempRealSpace = ifft(exp(-sum(dW)*1i*kSq).*currU);
                    % Full nonlinear step
                    currU = fft(exp(h*1i*absSq(tempRealSpace).^sigma).*tempRealSpace);
                case 8 % Fourier splitting scheme
                    % Full nonlinear step
                    tempRealSpace = ifft(currU);
                    temp = exp(h*1i*absSq(tempRealSpace).^sigma).*tempRealSpace;
                    % Full linear step
                    currU = exp(-sum(dW)*1i*kSq).*fft(temp);
                case 9 % Strang splitting scheme
                    % Half first step
                    temp = exp(-dW(1)*1i*kSq).*currU;
                    % Full nonlinear step
                    tempRealSpace = ifft(temp);
                    temp = fft(exp(h*1i*absSq(tempRealSpace).^sigma).*tempRealSpace);
                    % Half last step
                    currU = exp(-dW(2)*1i*kSq).*temp;
            end
            
            tempCurrU = ifft(currU);
            currL2Drift = abs(L2norm(tempCurrU) - u0L2Val);
            if currL2Drift > maxL2Drift(j)
                maxL2Drift(j) = currL2Drift;
            end
        end
    end
    maxL2DriftBatch(m,:) = maxL2Drift;
end

%% Vertical histogram plot
% Set y axis info and plot the information in form of histrograms
histWidh = 0.774;
schemePlotSpec = [0.1305,histWidh/numUsedSchemes,0.11,0.815];

logDiff = log10(maxL2DriftBatch);

yMax = ceil(max(logDiff(:)));
yMin = floor(min(logDiff(:)));

if yMin < -17
    yMin = -17;
end
if yMax > 2
    yMax = 2;
end

% Set the y ticks and labels
yAxisVector = yMin:0.2:yMax;

sideTitle = '$max_{n\in\{1,2,\ldots,N\}}\left(log(| ~ ||u_n||_{L^2}^2-||u_0||_{L^2}^2 ~ |)\right)$';
vertHistPlot(logDiff,schemes.shortNames(schemeIndexMat(:,2)),yAxisVector,sideTitle,schemePlotSpec)
set(gcf, 'Position', get(0, 'Screensize'));
%% Misc. functions
% Squared absolute value
% Also seen in absSq.m
function ret = absSq(u)
    ret = real(u.*conj(u));
end

% Crank-Nicolson non-linear function
% Also seen in CNNonLinG.m
function G = CNnonLin(u,v,sigma)
if sigma == 1
    G = (absSq(u)+absSq(v)).*(u+v)/(2*(sigma+1));
elseif sigma == 2
    G = (absSq(u).^2 +...
        absSq(u).*absSq(v) +...
        absSq(v).^2).*(u+v)/(2*(sigma+1));
elseif sigma == 3
    G = (absSq(u).^3 +...
        absSq(u).^2.*absSq(v) +...
        absSq(u).*absSq(v).^2+...
        absSq(v).^3).*(u+v)/(2*(sigma+1));
elseif sigma == 4
    G = (absSq(u).^4 +...
        absSq(u).^3.*absSq(v) +...
        absSq(u).^2.*absSq(v).^2+...
        absSq(u).*absSq(v).^3+...
        absSq(v).^4).*(u+v)/(2*(sigma+1));
else
    % This function is designed to avoid divison by zero
    G = (absSq(u).^(sigma+1)-absSq(v).^(sigma+1))./...
        (absSq(u)-absSq(v) + eps).*(u+v)/(2*(sigma+1));
end
end

% Vertical histogram plot function
% Also seen in vertHistPlot.m
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
    set(gca,'FontSize',25)

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

% Integral approximation using the trapezoidal method
% Also seen in trapezoidalIntegral.m
function [integral] = trapezoidalIntegral(dx,f)
% This function will approximate the integral int(f(x)) given the input:
%   dx  - Space step size
%   f   - A vector of f(x), where x is a vector such that x(i)-x(i-1)=dx
integral = dx/2*sum((f(1:(end-1)) + f(2:end)));
end