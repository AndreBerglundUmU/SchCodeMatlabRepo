function [ NStar, varargout ] = NStarSolverPS(currU,a,h,sigma)
crit = true;
K = length(currU);
NStar = currU;
i = 1;
while crit && i<120
    oldNStar = NStar;
    tempNStar = a+h/2*NStar;
    NStar = fft(1i*nonLin(ifft(tempNStar),sigma));
    crit = norm((oldNStar-NStar)./K,2) > eps;
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

