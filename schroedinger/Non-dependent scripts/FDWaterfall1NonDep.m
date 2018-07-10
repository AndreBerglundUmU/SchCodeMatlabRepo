%% Set initial info and function
initSeed = 1;
% Time and area
L = 15; XInt = [-L,L];
T = 1; TInt = [0,T];

% Numerical precision, number of points
N = 2^9; % Time
M = 2^7; % Space

h = (TInt(2)-TInt(1))/N;
dx = (XInt(2)-XInt(1))/M;

alpha = 1;
q = 1;
c = 1;
u0Fun = @(x) sqrt(2*alpha/q)*exp(0.5*1i*c*x).*sech(sqrt(alpha)*x);
% u0Fun = @(x) exp(-2*x.^2);

sigma = 1;
G = @(u) absSq(u).^sigma.*u;

% See finiteDifferenceSPMatrix.m
% finiteDifferenceSPMatrix(3,2,M+1)/dx^2
% Construct second order periodic derivative matrix
diagLine = ones(M+1,1);
FDMatSq = spdiags([-2*diagLine,diagLine,diagLine],[0,-1,1],M+1,M+1)/dx^2;
% Enforce Dirschlet boundary conditions
FDMatSq(1,:) = zeros(size(FDMatSq(1,:)));
FDMatSq(end,:) = zeros(size(FDMatSq(end,:)));
FDMatSq(:,1) = zeros(size(FDMatSq(:,1)));
FDMatSq(:,end) = zeros(size(FDMatSq(:,end)));

t = linspace(TInt(1),TInt(2),N+1);
x = XInt(1) + dx*(0:M)';
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
    currU = u0FunVal;
    currScheme = schemeIndexMat(j,2);
    queryIndex = 2;
    for i = 1:N
        dW = W(i,:);
        %% Scheme and query calculations
        switch currScheme
            case 1 % Forwards/explicit Euler scheme
                % Will perform one calculation in excess
                K = M+1;
                A = speye(K) + 1i*sum(dW)/2*FDMatSq;
                B = speye(K) - 1i*sum(dW)/2*FDMatSq;
                nextU = B\(A*currU + 1i*h*G(currU));
            case 2 % Backwards/implicit Euler scheme
                crit = true;
                K = M+1;
                nextU = currU;
                A = speye(K) + 1i*sum(dW)/2*FDMatSq;
                B = speye(K) - 1i*sum(dW)/2*FDMatSq;
                AcurrU = A*currU;
                
                noRounds = 1;
                while crit && noRounds<120
                    tempU = nextU;
                    nextU = B\(AcurrU + 1i*h*G(nextU));
                    crit = norm((tempU-nextU)./K,2) > eps;
                    noRounds = noRounds+1;
                end
                currU = nextU;
            case 3 % Midpoint scheme
                crit = true;
                K = M+1;
                nextU = currU;
                A = speye(K) + 1i*sum(dW)/2*FDMatSq;
                B = speye(K) - 1i*sum(dW)/2*FDMatSq;
                AcurrU = A*currU;
                
                noRounds = 1;
                while crit && noRounds<120
                    tempU = nextU;
                    nextU = B\(AcurrU + 1i*h*G((currU+nextU)/2));
                    crit = norm((tempU-nextU)./K,2) > eps;
                    noRounds = noRounds+1;
                end
                currU = nextU;
            case 4 % Crank-Nicolson scheme
                crit = true;
                K = M+1;
                nextU = currU;
                A = speye(K) + 1i*sum(dW)/2*FDMatSq;
                B = speye(K) - 1i*sum(dW)/2*FDMatSq;
                AcurrU = A*currU;
                
                noRounds = 1;
                while crit && noRounds<120
                    tempU = nextU;
                    nextU = B\(AcurrU + 1i*h*CNnonLin(currU,nextU,sigma));
                    crit = norm((tempU-nextU)./K,2) > eps;
                    noRounds = noRounds+1;
                end
                currU = nextU;
            case 5 % Explitic exponential scheme
                currU = expm(1i*sum(dW)*FDMatSq)*(currU + 1i*h*G(currU));
            case 6 % Symmetric exponential scheme
                % Implicit solving for NStar
                a = expm(1i*dW(1)*FDMatSq)*currU;
                crit = true;
                K = M+1;
                NStar = currU;
                noRounds = 1;
                while crit && noRounds<120
                    oldNStar = NStar;
                    tempNStar = a+h/2*NStar;
                    NStar = 1i*G(tempNStar);
                    crit = norm((oldNStar-NStar)./K,2) > eps;
                    noRounds = noRounds+1;
                end
                % Last step
                currU = expm(1i*sum(dW)*FDMatSq)*currU + h*expm(1i*dW(2)*FDMatSq)*NStar;
            case 7 % Lie-Trotter splitting scheme
                % Full linear step
                temp = expm(1i*sum(dW)*FDMatSq)*currU;
                % Full nonlinear step
                currU = exp(h*1i*absSq(temp).^sigma).*temp;
            case 8 % Fourier splitting scheme
                % Full nonlinear step
                temp = exp(h*1i*absSq(currU).^sigma).*currU;
                % Full linear step
                currU = expm(1i*sum(dW)*FDMatSq)*temp;
            case 9 % Strang splitting scheme
                % Half first step
                temp = expm(1i*dW(1)*FDMatSq)*currU;
                % Full nonlinear step
                temp = exp(h*1i*absSq(temp).^sigma).*temp;
                % Half last step
                currU = expm(1i*dW(2)*FDMatSq)*temp;
        end
        
        if mod(i,timeScalingFactor) == 0
            runStorage{j}(queryIndex,:) = currU(spaceVec);
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