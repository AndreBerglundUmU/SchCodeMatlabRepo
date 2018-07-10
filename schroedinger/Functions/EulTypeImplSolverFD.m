function [ nextU, varargout ] = EulTypeImplSolverFD(currU,FDMatSq,dW,h,G)
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