function [ nextU, varargout ] = EulTypeImplSolverFD(currU,FDMatSq,dW,h,G)
% EulTypeImplSolverFD  - Returns an implicitly solved u_(n+1) for Euler
%                        type schemes for the stochastic Schroedingers
%                        equation, with the finite difference approximation.
% Syntax: [ nextU, varargout ] = EulTypeImplSolverFD(currU,FDMatSq,dW,h,G)
%
% Input:
% currU - A vector containing u_n in Fourier space.
% a     - A vector containing exp(-dW(1)*1i*kSq).*currU.
% h     - The time step size.
% sigma - The scalar controlling the non-linearity.
%
% Output:
% nextU     - A vector of length M containing u_{n+1}.
% varargout - If asked for, a boolean value revealing whether the implicit
%             calculation converged in 120 fixed point iterations or not.
%
% Non-standard dependencies: None.
% See also: Any accompanying script for example usage.
%           makePSSchroedSchemes.m

crit = true;
K = length(currU);
nextU = currU;
A = speye(K) + 1i*dW/2*FDMatSq;
B = speye(K) - 1i*dW/2*FDMatSq;
AcurrU = A*currU;
i = 1;

while crit && i<120
    tempU = nextU;
    nextU = B\(AcurrU + 1i*h*G(currU,nextU));
    crit = norm((tempU-nextU)./K,2) > eps;
    i = i+1;
end
if nargout == 2
    if i == 120
        varargout{1} = true;
    else
        varargout{1} = false;
    end
end
end