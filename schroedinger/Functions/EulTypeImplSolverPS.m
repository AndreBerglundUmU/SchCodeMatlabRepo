function [ nextU, varargout ] = EulTypeImplSolverPS(currU,kSq,dW,h,G)
% EulTypeImplSolverPS   - Implicit solver for PS Euler type schemes.
%
% EulTypeImplSolverPS   - Returns an implicitly solved u_{n+1} for Euler
%                         type schemes for the stochastic Schroedingers
%                         equation, with the pseudospectral approximation.
% Syntax: [ nextU, varargout ] = EulTypeImplSolverPS(currU,kSq,dW,h,G)
%
% Input:
% currU     - A vector of length M containing u_n in Fourier space.
% kSq       - A vector of length M containing the squared Fourier modes
% dW        - A scalar value of the Brownian motion.
% h         - The time step size.
% G         - A function handle G(u_n,u_{n+1}) approximating the nonlinearity.
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
realSpaceCurrU = ifft(currU);
realSpaceNextU = realSpaceCurrU;

a = 1i*dW/2*kSq;
b = (1 - a)./(1 + a).*currU;
c = 1i*h./(1 + a);

i = 1;
while crit && i<120
    tempU = nextU;
    nextU = b + c.*fft(G(realSpaceCurrU,realSpaceNextU));
    realSpaceNextU=ifft(nextU);
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