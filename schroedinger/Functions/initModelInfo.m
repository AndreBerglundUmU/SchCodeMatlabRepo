function ret = initModelInfo(N,TInt,M,XInt,sigma,per)
    ret.N = N;
    ret.M = M;
    ret.h = (TInt(2)-TInt(1))/N;
    ret.dx = (XInt(2)-XInt(1))/M;

    % If periodic, PS. If not, FD
    if per
        ret.x = XInt(1) + ret.dx*(0:M-1);
        ret.k = 2*pi/(XInt(2)-XInt(1))*[0:M/2-1, 0, -M/2+1:-1];
        % Retrieve schemes
        ret.schemes = makePSSchroedSchemes(ret.k.^2,ret.h,sigma);
    else
        ret.x = XInt(1) + ret.dx*(0:M)';
        % Retrieve second order periodic derivative matrix
        FDMatSq = finiteDifferenceSpMatrix(3,2,M+1)/ret.dx^2;
        % Enforce Dirschlet boundary conditions
        FDMatSq(1,:) = zeros(size(FDMatSq(1,:)));
        FDMatSq(end,:) = zeros(size(FDMatSq(end,:)));
        FDMatSq(:,1) = zeros(size(FDMatSq(:,1)));
        FDMatSq(:,end) = zeros(size(FDMatSq(:,end)));
        ret.FDMatSq = FDMatSq;
        % Retrieve schemes
        ret.schemes = makeFDSchroedSchemes(ret.FDMatSq,ret.h,sigma);
    end
end