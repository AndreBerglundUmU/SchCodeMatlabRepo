addpath('Functions')
%% Set initial info and function
initSeed = 1;
% Time and area
XInt = [0,2*pi];
T = 1; TInt = [0,T];

% Numerical precision, number of points
N = 2^13; % Time
M = 2^13; % Space
sigma = 1;

u0Fun = @(x) 1./(2+sin(x).^2);
per = true;

modelInfo = initModelInfo(N,TInt,M,XInt,sigma,per);

u0FunVal = u0Fun(modelInfo.x);

numAvailableSchemes = length(modelInfo.schemes.fun); % 9 for this

schemesUsed = false(numAvailableSchemes,1);
% schemesUsed(1) = true; % FEul
% schemesUsed(2) = true; % BEul
% schemesUsed(3) = true; % MP
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
W = randn(N,2)*sqrt(modelInfo.h/2);
for j = 1:numUsedSchemes
    currU = fft(u0FunVal);
    currScheme = schemeIndexMat(j,2);
    queryIndex = 2;
    for i = 1:N
        dW = W(i,:);
        %% Scheme and query calculations
        currU = modelInfo.schemes.fun{currScheme}(currU,dW);
        
        tempCurrU = ifft(currU);
        if mod(i,timeScalingFactor) == 0
            runStorage{j}(queryIndex,:) = tempCurrU(spaceVec);
            queryIndex = queryIndex + 1;
        end
    end
end
%% Waterfall plots
for i = 1:numUsedSchemes
    figure
    waterFallPlot(runStorage{i},XInt,length(spaceVec),T,length(timeVec))
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'FontSize',35)
    pause(1)
%     printToEPS(gcf,['PSWaterfall' schemeShortNames{i}])
%     printToPDF(gcf,['PSWaterfall' schemes.shortNames{schemeIndexMat(i,2)}])
end