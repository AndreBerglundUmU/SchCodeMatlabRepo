function [ NStar, varargout ] = NStarSolverFD(currU,a,h,sigma)
% NStarSolverFD  -  Returns an implicitly solved N* for the symmetric
%                   exponential scheme for the stochastic Schroedingers
%                   equation, with the finite difference approximation.
% Syntax: [ NStar, varargout ] = NStarSolverPS(currU,a,h,sigma)
%
% Input:
% currU - A vector containing u_n in physical space.
% a     - A vector containing expm(1i*dW(1)*FDMatSq)*currU.
% h     - The time step size.
% sigma - The scalar controlling the non-linearity.
%
% Output:
% NStar     - A vector of length M containing NStar.
% varargout - If asked for, a boolean value revealing whether the implicit
%             calculation converged in 120 fixed point iterations or not.
%
% Non-standard dependencies: nonLin.m.
% See also: Any accompanying script for example usage.
%           makeFDSchroedSchemes.m

crit = true;
K = length(currU);
NStar = currU;
i = 1;
while crit && i<120
    oldNStar = NStar;
    tempNStar = a+h/2*NStar;
    NStar = 1i*nonLin(tempNStar,sigma);
    crit = norm((oldNStar-NStar)./K,2) > eps;
    i = i+1;
end
% Return convergence result if asked for
if nargout == 2
    if i == 120
        varargout{1} = true;
    else
        varargout{1} = false;
    end
end
end