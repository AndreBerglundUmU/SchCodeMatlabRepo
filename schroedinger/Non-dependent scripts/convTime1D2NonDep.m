%% Set initial info and function
initSeed = 1;
batchSize = 2000;
confLevel = 0.95;
% Time and area
L = 20; XInt = [-L,L];
T = 1; TInt = [0,T];

% Numerical precision, number of points
refN = 2^17;
refh = (TInt(2)-TInt(1))/refN;
N = 2.^(5:13); % Time
h = (TInt(2)-TInt(1))./N;
numN = length(N);
M = 2^11; % Space
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
L2DiffBatch = zeros(batchSize,length(N),numUsedSchemes);

%% Perform calculations
rng(initSeed,'twister')

parfor m = 1:batchSize
    % Load the broadcast variables to internal
    internalSchemeIndexMat = schemeIndexMat;
    internalN = N;
    internalh = h;
    
    % Calculate reference solution
    refW = randn(refN,2)*sqrt(refh/2);
    currU = fft(u0FunVal);
    for i = 1:refN
        dW = refW(i,:);
        % Using Strang splitting scheme as reference solution
        % Half first step
        temp = exp(-dW(1)*1i*kSq).*currU;
        % Full nonlinear step
        tempRealSpace = ifft(temp);
        temp = fft(exp(refh*1i*absSq(tempRealSpace).^sigma).*tempRealSpace);
        % Half last step
        currU = exp(-dW(2)*1i*kSq).*temp;
    end
    refSol = ifft(currU);
    
    tempBatch = zeros(length(N),numUsedSchemes);
    for n = 1:length(N)
        % Retrieve information regarding current number of steps
        currN = internalN(n);
        currh = internalh(n);
        scalingFactor = refN/currN;
        
        l2Difference = zeros(1,numUsedSchemes);
        for j = 1:numUsedSchemes
            currU = fft(u0FunVal);
            currScheme = internalSchemeIndexMat(j,2);
            currIndex = 0;
            dW = zeros(2,1);
            for i = 1:currN
                % Calculate the coarser Brownian motion increments
                indexList = (currIndex+1):(currIndex+scalingFactor/2);
                dW(1) = sum(sum(refW(indexList,:)));
                indexList = (currIndex+scalingFactor/2+1):(currIndex+scalingFactor);
                dW(2) = sum(sum(refW(indexList,:)));
                currIndex = currIndex + scalingFactor;
                %% Scheme and query calculations
                % See makePSSchroedSchemes.m for schemes
                switch currScheme
                    case 1 % Forwards/explicit Euler scheme
                        realSpaceCurrU = ifft(currU);
                        a = 1i*sum(dW)/2*kSq;
                        b = (1 - a)./(1 + a).*currU;
                        c = 1i*currh./(1 + a);
                        currU = b + c.*fft(G(realSpaceCurrU));

                    case 2 % Backwards/implicit Euler scheme
                        % Implicit solving for next currU
                        crit = true;
                        nextU = currU;
                        realSpaceCurrU = ifft(currU);
                        realSpaceNextU = realSpaceCurrU;
                        
                        a = 1i*sum(dW)/2*kSq;
                        b = (1 - a)./(1 + a).*currU;
                        c = 1i*currh./(1 + a);
                        
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
                        c = 1i*currh./(1 + a);
                        
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
                        c = 1i*currh./(1 + a);
                        
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
                        currU = exp(-sum(dW)*1i*kSq).*(currU + 1i*currh*fft(G(ifft(currU))));
                    case 6 % Symmetric exponential scheme
                        % Implicit solving for NStar
                        a = exp(-dW(1)*1i*kSq).*currU;
                        crit = true;
                        NStar = currU;
                        noRounds = 1;
                        while crit && noRounds<120
                            oldNStar = NStar;
                            tempNStar = a+currh/2*NStar;
                            NStar = fft(1i*G(ifft(tempNStar)));
                            crit = norm((oldNStar-NStar)./M,2) > eps;
                            noRounds = noRounds+1;
                        end
                        % Last step
                        currU = exp(-sum(dW)*1i*kSq).*currU + currh*exp(-dW(2)*1i*kSq).*NStar;
                    case 7 % Lie-Trotter splitting scheme
                        % Full linear step
                        tempRealSpace = ifft(exp(-sum(dW)*1i*kSq).*currU);
                        % Full nonlinear step
                        currU = fft(exp(currh*1i*absSq(tempRealSpace).^sigma).*tempRealSpace);
                    case 8 % Fourier splitting scheme
                        % Full nonlinear step
                        tempRealSpace = ifft(currU);
                        temp = exp(currh*1i*absSq(tempRealSpace).^sigma).*tempRealSpace;
                        % Full linear step
                        currU = exp(-sum(dW)*1i*kSq).*fft(temp);
                    case 9 % Strang splitting scheme
                        % Half first step
                        temp = exp(-dW(1)*1i*kSq).*currU;
                        % Full nonlinear step
                        tempRealSpace = ifft(temp);
                        temp = fft(exp(currh*1i*absSq(tempRealSpace).^sigma).*tempRealSpace);
                        % Half last step
                        currU = exp(-dW(2)*1i*kSq).*temp;
                end
            end
            l2Difference(j) = L2norm(refSol-ifft(currU));
        end
        tempBatch(n,:) = l2Difference;
    end
    L2DiffBatch(m,:,:) = tempBatch;
    m
end

%% Plot convergence

% Add support lines
numSupportLines = 2;
translConst = [4 3];
    
logInArg = cell(0);
legendInArg = cell(numSupportLines,1);
for i = 1:numSupportLines
    logInArg(end+1:end+3) = {N,translConst(i)*N.^-i,'--k'};
    legendInArg{i} = sprintf('N^{-%d}',i);
end

figure
hold on
for i = 1:numUsedSchemes
    tempMean = mean(L2DiffBatch(:,:,i));
    
    % Calculate CI based on normal distribution
    %         upperBar = abs(tempMean-prctile(l2Storage(i,:,:),confLevel*100,3));
    %         lowerBar = abs(tempMean-prctile(l2Storage(i,:,:),(1-confLevel)*100,3));
    %         errorbar(N,tempMean,lowerBar,upperBar,'-x');
    perc = abs(norminv((1-confLevel)/2,0,1));
    tempStd = std(L2DiffBatch(:,:,i));
    switch schemeIndexMat(i,2)
        case 1
            line = '-^b';
        case 2
            line = '-vm';
        case 5
            line = '-pc';
        case 6
            line = '-og';
        case 10
            line = '-or';
        otherwise
            line = '-x';
    end
    errorbar(N,tempMean,perc*tempStd/sqrt(batchSize),line,'LineWidth',1.5);
end

set(gca,'xscale','log');
set(gca,'yscale','log');
loglog(logInArg{:})
legend(refSchemes.shortNames{schemesUsed},legendInArg{:})
hold off

xlabel('N')
ylabel('$E[||u_N - u^*(T)||_{L^2}^2]$','Interpreter','latex')
set(gca,'FontSize',25)
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

function [integral] = trapezoidalIntegral(dx,f)
% This function will approximate the integral int(f(x)) given the input:
%   dx  - Space step size
%   f   - A vector of f(x), where x is a vector such that x(i)-x(i-1)=dx
integral = dx/2*sum((f(1:(end-1)) + f(2:end)));
end