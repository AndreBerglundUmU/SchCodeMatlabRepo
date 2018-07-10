%% Set initial info and function
initSeed = 1;
% Time and area
% L = 8; XInt = [-L,L];
XInt = [0,2*pi];
T = 1; TInt = [0,T];

% Numerical precision, number of points
N = 2^13; % Time
M = 2^13; % Space

h = (TInt(2)-TInt(1))/N;
dx = (XInt(2)-XInt(1))/M;

alpha = 1;
q = 1;
c = 1;
u0Fun = @(x) sqrt(2*alpha/q)*exp(0.5*1i*c*x).*sech(sqrt(alpha)*x);
% u0Fun = @(x) exp(-2*x.^2);

sigma = 1;
G = @(u) absSq(u).^sigma.*u;

k = 2*pi/(XInt(2)-XInt(1))*[0:M/2-1, 0, -M/2+1:-1];
kSq = k.^2;

t = linspace(TInt(1),TInt(2),N+1);
x = XInt(1) + dx*(0:M-1);
u0FunVal = u0Fun(x);

% See makeFDSchroedSchemes.m for schemes
numAvailableSchemes = 9; % 9 for this

schemesUsed = false(numAvailableSchemes,1);
% schemesUsed(1) = true; % FEul
% schemesUsed(2) = true; % BEul
% schemesUsed(3) = true; % MEul
schemesUsed(4) = true; % CN
schemesUsed(5) = true; % EExp
% schemesUsed(6) = true; % SExp
schemesUsed(7) = true; % LTSpl
% schemesUsed(8) = true; % FSpl
% schemesUsed(9) = true; % SSpl

numUsedSchemes = sum(schemesUsed);
schemeIndexMat = [(1:numUsedSchemes)' , find(schemesUsed)];

%% Query storage
% Query fidelity
timeFidelity = 2^7;
spaceFidelity = 2^9;

if N > timeFidelity
    timeScalingFactor = N / timeFidelity;
    timeVec = 0:timeScalingFactor:N;
else
    timeScalingFactor = 1;
    timeVec = 0:N;
end
if M > spaceFidelity
    spaceScalingFactor = M / spaceFidelity;
    spaceVec = 1:spaceScalingFactor:M;
else
    spaceVec = 1:M+1;
end

runStorage = cell(numUsedSchemes,1);
for i = 1:numUsedSchemes
    runStorage{i} = zeros(length(timeVec),length(spaceVec));
    runStorage{i}(1,:) = u0FunVal(spaceVec);
end
%% Perform calculations
rng(initSeed,'twister')
W = randn(N,2)*sqrt(h/2);
for j = 1:numUsedSchemes
    currU = fft(u0FunVal);
    currScheme = schemeIndexMat(j,2);
    queryIndex = 2;
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
        if mod(i,timeScalingFactor) == 0
            runStorage{j}(queryIndex,:) = tempCurrU(spaceVec);
            queryIndex = queryIndex + 1;
        end
    end
end
%% Waterfall plots
for i = 1:numUsedSchemes
    % Also seen in waterFallPlot.m
    figure
    xVec = linspace(XInt(1),XInt(2),length(spaceVec));
    tVec = linspace(0,T,length(timeVec));
    temp = waterfall(xVec,tVec,absSq(runStorage{i}));
    set(temp,'LineWidth',1.5)
    ylabel('t','FontSize',20);
    xlabel('x','FontSize',20);
    zlabel('$|u_n|^2$','Rotation',0,'HorizontalAlignment','right','Interpreter','latex','FontSize',20);
    axis tight
    
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'FontSize',25)
    pause(1)
end
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