function ret = nonLin(u,sigma)
% A short function for the nonlinearity in the given Schroedinger equation.
% ret =  nonLin(u,sigma) = absSq(u).^sigma.*u;
    ret = absSq(u).^sigma.*u;
end