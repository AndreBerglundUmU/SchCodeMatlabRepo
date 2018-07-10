function [ nextU, varargout ] = EulTypeImplSolverPS(currU,kSq,dW,h,G)
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