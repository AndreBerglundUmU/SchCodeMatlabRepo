function ret = L2norm(u,dx,per)
    if per
        ret = trapezoidalIntegral(dx,absSq(u)) + dx*(absSq(u(end))+absSq(u(1)))/2;
    else
        ret = trapezoidalIntegral(dx,absSq(u));
    end
end