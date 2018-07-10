function schemes = makePSSchroedSchemes(kSq,h,sigma)
    schemes.fun = cell(9,1);
        schemes.fun{1} = @(currU,dW) FEul(currU,kSq,dW,h,sigma);
        schemes.fun{2} = @(currU,dW) BEul(currU,kSq,dW,h,sigma);
        schemes.fun{3} = @(currU,dW) MEul(currU,kSq,dW,h,sigma);
        schemes.fun{4} = @(currU,dW) CN(currU,kSq,dW,h,sigma);
        schemes.fun{5} = @(currU,dW) EExp(currU,kSq,dW,h,sigma);
        schemes.fun{6} = @(currU,dW) SExp(currU,kSq,dW,h,sigma);
        schemes.fun{7} = @(currU,dW) LTSpl(currU,kSq,dW,h,sigma);
        schemes.fun{8} = @(currU,dW) FSpl(currU,kSq,dW,h,sigma);
        schemes.fun{9} = @(currU,dW) SSpl(currU,kSq,dW,h,sigma);
    schemes.longNames = cell(9,1);
        schemes.longNames{1} = 'PS explicit "Euler"';
        schemes.longNames{2} = 'PS implicit "Euler"';
        schemes.longNames{3} = 'PS midpoint';
        schemes.longNames{4} = 'PS Crank-Nicolson';
        schemes.longNames{5} = 'PS exlicit exponential';
        schemes.longNames{6} = 'PS symmetric exponential';
        schemes.longNames{7} = 'PS Lie splitting';
        schemes.longNames{8} = 'PS Fourier splitting';
        schemes.longNames{9} = 'PS Strang splitting';
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

function nextU = FEul(currU,kSq,dW,h,sigma)
    % Will perform one calculation in excess
    % nextU = EulTypeImplSolverPS(currU,kSq,sum(dW),h,@(un,unext) G(un,sigma));
    
    realSpaceCurrU = ifft(currU);
    a = 1i*sum(dW)/2*kSq;
    b = (1 - a)./(1 + a).*currU;
    c = 1i*h./(1 + a);
    nextU = b + c.*fft(nonLin(realSpaceCurrU,sigma));
end

function nextU = BEul(currU,kSq,dW,h,sigma)
    nextU = EulTypeImplSolverPS(currU,kSq,sum(dW),h,@(un,unext) nonLin(unext,sigma));
end

function nextU = MEul(currU,kSq,dW,h,sigma)
    nextU = EulTypeImplSolverPS(currU,kSq,sum(dW),h,@(un,unext) nonLin((un+unext)/2,sigma));
end

function nextU = CN(currU,kSq,dW,h,sigma)
    nextU = EulTypeImplSolverPS(currU,kSq,sum(dW),h,@(un,unext) CNnonLin(un,unext,sigma));
end

function nextU = EExp(currU,kSq,dW,h,sigma)
    nextU = exp(-sum(dW)*1i*kSq).*(currU + 1i*h*fft(nonLin(ifft(currU),sigma)));
end

function nextU = SExp(currU,kSq,dW,h,sigma)
    a = exp(-dW(1)*1i*kSq).*currU;
    NStar = NStarSolverPS(currU,a,h,sigma);
    nextU = exp(-sum(dW)*1i*kSq).*currU + h*exp(-dW(2)*1i*kSq).*NStar;
end

function nextU = LTSpl(currU,kSq,dW,h,sigma)
    % Full linear step
    tempRealSpace = ifft(exp(-sum(dW)*1i*kSq).*currU);
    % Full nonlinear step
    nextU = fft(exp(h*1i*absSq(tempRealSpace).^sigma).*tempRealSpace);
end

function nextU = FSpl(currU,kSq,dW,h,sigma)
    % Full nonlinear step
    tempRealSpace = ifft(currU);
    temp = exp(h*1i*absSq(tempRealSpace).^sigma).*tempRealSpace;
    % Full linear step
    nextU = exp(-sum(dW)*1i*kSq).*fft(temp);
end

function nextU = SSpl(currU,kSq,dW,h,sigma)
    % Half first step
    temp = exp(-dW(1)*1i*kSq).*currU;
    % Full nonlinear step
    tempRealSpace = ifft(temp);
    temp = fft(exp(h*1i*absSq(tempRealSpace).^sigma).*tempRealSpace);
    % Half last step
    nextU = exp(-dW(2)*1i*kSq).*temp;
end