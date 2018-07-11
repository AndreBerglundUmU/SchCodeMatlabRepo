addpath('Functions')
%% Set initial info and function
initSeed = 1;
% Time and area
L = 30; XInt = [-L,L];
T = 1; TInt = [0,T];

% Numerical precision, number of points
N = 2^9; % Time
M = 2^8; % Space
sigma = 1;

u0Fun = @(x) exp(-2*x.^2);
per = false;

modelInfo = initModelInfo(N,TInt,M,XInt,sigma,per);

u0FunVal = u0Fun(modelInfo.x);

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
% Query fidelity
timeFidelity = 2^5;
if N > timeFidelity
    scalingFactor = N / timeFidelity;
    timeVec = 0:scalingFactor:N;
else
    scalingFactor = 1;
    timeVec = 0:N;
end

normStorage = cell(numUsedSchemes,1);
for i = 1:numUsedSchemes
    normStorage{i} = zeros(length(timeVec),1);
    normStorage{i}(1) = L2norm(u0FunVal,modelInfo.dx,per);
end

%% Perform calculations
rng(initSeed,'twister')
W = randn(N,2)*sqrt(modelInfo.h/2);
for j = 1:numUsedSchemes
    currU = u0FunVal;
    currScheme = schemeIndexMat(j,2);
    queryIndex = 2;
    for i = 1:N
        dW = W(i,:);
        %% Scheme and query calculations
        currU = modelInfo.schemes.fun{currScheme}(currU,dW);
        
        if mod(i,scalingFactor) == 0
            normStorage{j}(queryIndex) = L2norm(currU,modelInfo.dx,per);
            queryIndex = queryIndex + 1;
        end
    end
end

%% L2 norm evolution plot
inArg = cell(3*numUsedSchemes,1);
legendInArg = cell(numUsedSchemes,1);

coarseTVec = linspace(TInt(1),TInt(2),timeFidelity+1);

maxUpwardDrift = 0;
maxDownwardDrift = 0;
for i = 1:numUsedSchemes
    % Retrieve info of which scheme was ran
    usedScheme = schemeIndexMat(i,2);
    switch usedScheme
        case 1
            inArg(3*i-2:3*i) =  {coarseTVec,normStorage{i},'-^r'};
        case 2
            inArg(3*i-2:3*i) =  {coarseTVec,normStorage{i},'-vm'};
        case 5
            inArg(3*i-2:3*i) =  {coarseTVec,normStorage{i},'-db'};
        otherwise
            inArg(3*i-2:3*i) =  {coarseTVec,normStorage{i},'k'};
    end
    % Retrieve legend info
    legendInArg(i) = modelInfo.schemes.shortNames(usedScheme);
    % Retrieve extreme values
    driftNormVector = normStorage{i}-normStorage{i}(1);
    driftUpwardIndex = driftNormVector >= 0;
    
    schemeDriftUpward = max(driftNormVector(driftUpwardIndex));
    schemeDriftDownward = min(driftNormVector(~driftUpwardIndex));
    
    if maxUpwardDrift < schemeDriftUpward
        maxUpwardDrift = schemeDriftUpward;
    end
    if maxDownwardDrift > schemeDriftDownward
        maxDownwardDrift = schemeDriftDownward;
    end
end

% Make sure that the plot window won't touch the lines
if -maxDownwardDrift < 0.1*maxUpwardDrift
    maxDownwardDrift = -0.1*maxUpwardDrift;
elseif -0.1*maxDownwardDrift > maxUpwardDrift
    maxUpwardDrift = -0.1*maxDownwardDrift;
end

figure
plot(inArg{:},'LineWidth',1.5,'MarkerSize',13)
legend(legendInArg{:},'Location','northwest')
ylabel('$||u_n||_{L^2}^2$','Interpreter','latex')
xlabel('t');
set(gca,'FontSize',35)

axis([TInt normStorage{1}(1)+1.1*maxDownwardDrift normStorage{1}(1)+1.1*maxUpwardDrift])
set(gcf, 'Position', get(0, 'Screensize'));

pause(1)
% printToPDF(gcf,'FDEvol')