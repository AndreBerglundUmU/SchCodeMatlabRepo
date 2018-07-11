function ret = initModelInfo(N,TInt,M,XInt,sigma,per)
% initModelInfo  -  Returns a struct containing the 
% Syntax: ret = initModelInfo(N,TInt,M,XInt,sigma,per)
%
% Input:
% N     - The number of time points.
% TInt  - A vector of length 2 containing the time interval end points.
% M     - The number of spacial points.
% XInt  - A vector of length 2 containing the space interval end point
% sigma - The scalar controlling the non-linearity.
% per   - A boolean value declaring whether the problem is periodic or not.
%
% Output:
% ret   - A struct containing the fields
%         N - The number of time points.
%         M - The number of spacial points.
%         h - The time step size.
%         dx - The spacial step size.
%         x - The space vector
%         k - Only available if per. The Fourier mode vector.
%         FDMatSq - Only available if ~per. The second order differential.
%         schemes - A struct containing either the FD or PS schemes.
%
% Non-standard dependencies: makePSSchroedSchemes.m ,
%                            makeFDSchroedSchemes.m,
%                            finiteDifferenceSpMatrix.m.
% See also: Any accompanying script for example usage.
%           makePSSchroedSchemes.m
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