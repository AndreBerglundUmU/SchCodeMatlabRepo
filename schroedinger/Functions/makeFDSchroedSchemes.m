function schemes = makeFDSchroedSchemes(FDMatSq,h,sigma)
% makeFDSchroedSchemes  - A function which will construct the schemes for
%                         the stochastic Schroedingers equation with white
%                         noise dispersion, with the finite difference 
%                         approximation.
% Syntax: schemes = makeFDSchroedSchemes(FDMatSq,h,sigma)
%
% Input:
% FDMatSq   - A MxM matrix containing the finite difference approximation
%             of the second order derivative.
% h         - The time step size.
% sigma     - The scalar controlling the non-linearity.
%
% Output:
% schemes   - A struct containing the fields .fun, .longNames and
%             .shortNames. Currently 9 schemes are implemented.
%
% Non-standard dependencies: absSq.m, nonLin.m, CNnonLin.m,
%                            EulTypeImplSolverFD.m, NStarSolverFD.m.
% See also: Any accompanying script for example usage.

    schemes.fun = cell(9,1);
        schemes.fun{1} = @(currU,dW) FEul(currU,FDMatSq,dW,h,sigma);
        schemes.fun{2} = @(currU,dW) BEul(currU,FDMatSq,dW,h,sigma);
        schemes.fun{3} = @(currU,dW) MEul(currU,FDMatSq,dW,h,sigma);
        schemes.fun{4} = @(currU,dW) CN(currU,FDMatSq,dW,h,sigma);
        schemes.fun{5} = @(currU,dW) EExp(currU,FDMatSq,dW,h,sigma);
        schemes.fun{6} = @(currU,dW) SExp(currU,FDMatSq,dW,h,sigma);
        schemes.fun{7} = @(currU,dW) LTSpl(currU,FDMatSq,dW,h,sigma);
        schemes.fun{8} = @(currU,dW) FSpl(currU,FDMatSq,dW,h,sigma);
        schemes.fun{9} = @(currU,dW) SSpl(currU,FDMatSq,dW,h,sigma);
    schemes.longNames = cell(9,1);
        schemes.longNames{1} = 'FD explicit "Euler"';
        schemes.longNames{2} = 'FD implicit "Euler"';
        schemes.longNames{3} = 'FD midpoint';
        schemes.longNames{4} = 'FD Crank-Nicolson';
        schemes.longNames{5} = 'FD exlicit exponential';
        schemes.longNames{6} = 'FD symmetric exponential';
        schemes.longNames{7} = 'FD Lie splitting';
        schemes.longNames{8} = 'FD Fourier splitting';
        schemes.longNames{9} = 'FD Strang splitting';
    schemes.shortNames = cell(9,1);
        schemes.shortNames{1} = 'FEul';
        schemes.shortNames{2} = 'BEul';
        schemes.shortNames{3} = 'MP';
        schemes.shortNames{4} = 'CN';
        schemes.shortNames{5} = 'EExp';
        schemes.shortNames{6} = 'SExp';
        schemes.shortNames{7} = 'LTSpl';
        schemes.shortNames{8} = 'FSpl';
        schemes.shortNames{9} = 'SSpl';
end

function nextU = FEul(currU,FDMatSq,dW,h,sigma)
    % Using the standard form will perform one calculation in excess:
    % nextU = FDEulTypeImplSolver(currU,FDMatSq,sum(dW),h,@(un,unext) G(un,sigma));
    
    K = size(FDMatSq,1);
    A = speye(K) + 1i*sum(dW)/2*FDMatSq;
    B = speye(K) - 1i*sum(dW)/2*FDMatSq;
    nextU = B\(A*currU + 1i*h*nonLin(currU,sigma));
end

function nextU = BEul(currU,FDMatSq,dW,h,sigma)
    nextU = EulTypeImplSolverFD(currU,FDMatSq,sum(dW),h,@(un,unext) nonLin(unext,sigma));
end

function nextU = MEul(currU,FDMatSq,dW,h,sigma)
    nextU = EulTypeImplSolverFD(currU,FDMatSq,sum(dW),h,@(un,unext) nonLin((un+unext)/2,sigma));
end

function nextU = CN(currU,FDMatSq,dW,h,sigma)
    nextU = EulTypeImplSolverFD(currU,FDMatSq,sum(dW),h,@(un,unext) CNnonLin(un,unext,sigma));
end

function nextU = EExp(currU,FDMatSq,dW,h,sigma)
    nextU = expm(1i*sum(dW)*FDMatSq)*(currU + 1i*h*nonLin(currU,sigma));
end

function nextU = SExp(currU,FDMatSq,dW,h,sigma)
    a = expm(1i*dW(1)*FDMatSq)*currU;
    NStar = NStarSolverFD(currU,a,h,sigma);
    nextU = expm(1i*sum(dW)*FDMatSq)*currU + h*expm(1i*dW(2)*FDMatSq)*NStar;
end

function nextU = LTSpl(currU,FDMatSq,dW,h,sigma)
    % Full linear step
    temp = expm(1i*sum(dW)*FDMatSq)*currU;
    % Full nonlinear step
    nextU = exp(h*1i*absSq(temp).^sigma).*temp;
end

function nextU = FSpl(currU,FDMatSq,dW,h,sigma)
    % Full nonlinear step
    temp = exp(h*1i*absSq(currU).^sigma).*currU;
    % Full linear step
    nextU = expm(1i*sum(dW)*FDMatSq)*temp;
end

function nextU = SSpl(currU,FDMatSq,dW,h,sigma)
    % Half first step
    temp = expm(1i*dW(1)*FDMatSq)*currU;
    % Full nonlinear step
    temp = exp(h*1i*absSq(temp).^sigma).*temp;
    % Half last step
    nextU = expm(1i*dW(2)*FDMatSq)*temp;
end