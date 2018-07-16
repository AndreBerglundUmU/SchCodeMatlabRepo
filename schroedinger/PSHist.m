addpath('Functions')
%% Set initial info and function
initSeed = 1;
batchSize = 1000;
% Time and area
XInt = [0,2*pi];
T = 1; TInt = [0,T];

% Numerical precision, number of points
N = 2^10; % Time
M = 2^10; % Space
sigma = 1;

u0Fun = @(x) 1./(2+sin(x).^2);
per = true;

modelInfo = initModelInfo(N,TInt,M,XInt,sigma,per);
% Needed for parfor loop
dx = modelInfo.dx;
h = modelInfo.h;

u0FunVal = u0Fun(modelInfo.x);
u0L2Val = L2norm(u0FunVal,dx,per);

numAvailableSchemes = length(modelInfo.schemes.fun); % 9 for this

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
maxL2DriftBatch = zeros(batchSize,numUsedSchemes);

%% Perform calculations
rng(initSeed,'twister')
parfor m = 1:batchSize
    % Load the broadcast variables to internal
    internalSchemes = modelInfo.schemes;
    internalSchemeIndexMat = schemeIndexMat;
    
    W = randn(N,2)*sqrt(h/2);
    
    maxL2Drift = zeros(1,numUsedSchemes);
    for j = 1:numUsedSchemes
        currU = fft(u0FunVal);
        currScheme = internalSchemeIndexMat(j,2);
        for i = 1:N
            dW = W(i,:);
            %% Scheme and query calculations
            currU = internalSchemes.fun{currScheme}(currU,dW);
            
            tempCurrU = ifft(currU);
            currL2Drift = abs(L2norm(tempCurrU,dx,per) - u0L2Val);
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
vertHistPlot(logDiff,modelInfo.schemes.shortNames(schemeIndexMat(:,2)),yAxisVector,sideTitle,schemePlotSpec)
set(gcf, 'Position', get(0, 'Screensize'));
pause(1)
% printToPDF(gcf,'PSHist')