function ret = nonLin(u,sigma)
    ret = absSq(u).^sigma.*u;
end