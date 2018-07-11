addpath('Functions')
%% Set initial info and function
initSeed = 1;
batchSize = 2000;
confLevel = 0.95;
% Time and area
XInt = [0,2*pi];
T = 1; TInt = [0,T];

% Numerical precision, number of points
refN = 2^17;
refh = (TInt(2)-TInt(1))/refN;
N = 2.^(5:13); % Time
h = (TInt(2)-TInt(1))./N;
numN = length(N);
M = 2^11; % Space
dx = (XInt(2)-XInt(1))/M;

u0Fun = @(x) 1./(2+sin(x).^2);
per = true;

sigma = 1;

k = 2*pi/(XInt(2)-XInt(1))*[0:M/2-1, 0, -M/2+1:-1];
kSq = k.^2;

x = XInt(1) + dx*(0:M-1);
u0FunVal = u0Fun(x);
u0L2Val = L2norm(u0FunVal,dx,per);

refSchemes = makePSSchroedSchemes(kSq,refh,sigma);
numAvailableSchemes = length(refSchemes.fun); % 9 for this

schemesUsed = false(numAvailableSchemes,1);
% schemesUsed(1) = true; % FEul
% schemesUsed(2) = true; % BEul
schemesUsed(3) = true; % MP
schemesUsed(4) = true; % CN
schemesUsed(5) = true; % EExp
schemesUsed(6) = true; % SExp
schemesUsed(7) = true; % LTSpl
% schemesUsed(8) = true; % FSpl
schemesUsed(9) = true; % SSpl

numUsedSchemes = sum(schemesUsed);
schemeIndexMat = [(1:numUsedSchemes)' , find(schemesUsed)];
%% Query storage
L2DiffBatch = zeros(batchSize,length(N),numUsedSchemes);

%% Perform calculations
rng(initSeed,'twister')
parfor m = 1:batchSize
    % Load the broadcast variables to internal
    internalRefSchemes = refSchemes;
    internalSchemeIndexMat = schemeIndexMat;
    internalN = N;
    internalh = h;
    
    % Calculate reference solution
    refW = randn(refN,2)*sqrt(refh/2);
    currU = fft(u0FunVal);
    for i = 1:refN
        dW = refW(i,:);
        % Using Strang splitting scheme as reference solution
        currU = internalRefSchemes.fun{9}(currU,dW);
    end
    refSol = ifft(currU);
    
    tempBatch = zeros(length(N),numUsedSchemes);
    for n = 1:length(N)
        % Retrieve information regarding current number of steps
        currN = internalN(n);
        currh = internalh(n);
        scalingFactor = refN/currN;
        schemes = makePSSchroedSchemes(kSq,currh,sigma);
        
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
                currU = schemes.fun{currScheme}(currU,dW);
            end
            l2Difference(j) = L2norm(refSol-ifft(currU),dx,per);
        end
        tempBatch(n,:) = l2Difference;
    end
    L2DiffBatch(m,:,:) = tempBatch
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
set(gca,'FontSize',35)
set(gcf, 'Position', get(0, 'Screensize'));

pause(1)
% printToPDF(gcf,'convTime1d')