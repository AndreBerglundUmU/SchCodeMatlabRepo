%% Set initial info and function
initSeed = 4;
batchSize = 4;
% Time and area
L = 20*pi; XInt = [-L,L];
% XInt = [0,2*pi];
T = 5*10^(-5); TInt = [0,T];

% Numerical precision, number of points
N = 2^18; % Time
M = 2^21; % Space
% N = 2^29; % Time
% M = 2^21; % Space
per = true;
sigma = 4;

modelInfo = initModelInfo(N,TInt,M,XInt,sigma,per);
%%
% u0Fun = @(x) 1./(2+sin(x).^2);
u0Fun = @(x) 4.5*exp(-4*x.^2);

L2norm = @(currU,dx) trapezoidalIntegral(dx,absSq([currU currU(1)]));
H1norm = @(currU,k,dx) L2norm(ifft(currU),dx) + L2norm(ifft(1i*k.*currU),dx);

u0FunVal = u0Fun(modelInfo.x);
u0H1Val = H1norm(fft(u0FunVal),modelInfo.k,modelInfo.dx);

numAvailableSchemes = length(modelInfo.schemes.fun); % 9 for this

schemesUsed = false(numAvailableSchemes,1);
% schemesUsed(1) = true; % FEul
% schemesUsed(2) = true; % BEul
% schemesUsed(3) = true; % MEul
% schemesUsed(4) = true; % CN
% schemesUsed(5) = true; % EExp
% schemesUsed(6) = true; % SExp
schemesUsed(7) = true; % LTSpl
% schemesUsed(8) = true; % FSpl
% schemesUsed(9) = true; % SSpl

numUsedSchemes = sum(schemesUsed);
schemeIndexMat = [(1:numUsedSchemes)' , find(schemesUsed)];

%% Query storage
timeInfoBatch = zeros(batchSize,3,N+1,numUsedSchemes);
timeInfoBatch(:,3,1,:) = u0H1Val*ones(batchSize,1,numUsedSchemes);

%% Perform calculations
rng(initSeed,'twister')
test = tic;
parfor m = 1:batchSize
    timeInfo = zeros(3,N+1,numUsedSchemes);
    internalSchemeIndexMat = schemeIndexMat;
    internalModelInfo = modelInfo;
    for j = 1:numUsedSchemes
        currU = fft(u0FunVal);
        currScheme = internalSchemeIndexMat(j,2);
        t = 0;
        W = 0;
        latestH1Val = u0H1Val;

        for i = 1:N
            dW = randn(2,1)*sqrt(internalModelInfo.h);
            %% Scheme and query calculations
            %%
%             test = tic;
%             for m = 1:10000
%             % Full linear step
%             tempRealSpace = ifft(exp(-sum(dW)*1i*kSq).*currU);
%             % Full nonlinear step
%             nextU = fft(exp(h*1i*absSq(tempRealSpace).^(sigma)).*tempRealSpace);
%             end
%             t1 = toc(test)

%             test = tic;
%             for m = 1:1000
            currU = internalModelInfo.schemes.fun{currScheme}(currU,dW);
%             end
%             t2 = toc(test)
%             t1/t2
            %%
            H1Val = H1norm(currU,internalModelInfo.k,internalModelInfo.dx);
            % Check if refinement is necessary
%             if (H1Val > 1.1*latestH1Val) && N < 2^26
%                 currTAndH1(:,i) = [t,H1Val];
%                 i = i + 1;
% 
%                 N = 2*N;
%                 M = 2*M;
%                 modelInfo = initModelInfo(N,TInt,M,XInt,sigma,per);
% 
%                 currU = periodicRefineVector(currU);
%             end

            t = t + modelInfo.h;
            W = W + sum(dW);

            timeInfo(1,i+1,j) = t;
            timeInfo(2,i+1,j) = W;
            timeInfo(3,i+1,j) = H1Val;

            disp([m t/T])
        end
    end
    timeInfoBatch(m,:,:,:) = timeInfo;
end
%%
figure
plot(squeeze(timeInfoBatch(1,1,:,1)),squeeze(timeInfoBatch(:,3,:,:)))
set(gca,'yscale','log');
xlabel('t')
ylabel('$||u_N||_{H^1}^2$','Interpreter','latex')
set(gca,'FontSize',35)
set(gcf, 'Position', get(0, 'Screensize'));

pause(1)
% printToPDF(gcf,'H1Explosion')