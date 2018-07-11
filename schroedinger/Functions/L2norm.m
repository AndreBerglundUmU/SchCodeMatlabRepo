function ret = L2norm(u,dx,per)
% L2norm - An approximation of the L2 norm, or the squared absolute integral
% Syntax: ret = L2norm(u,dx,per)
%
% Input:
% u     - A vector containing the function u in physical space
% dx    - The space step size.
% per   - A boolean value declaring whether u periodic or not.
%
% Output:
% ret - int abs(u)^2 dx
%
% Non-standard dependencies: absSq.m.
% See also: PSHist.m for example usage.
%           H1norm.m
    if per
        ret = trapezoidalIntegral(dx,absSq(u)) + dx*(absSq(u(end))+absSq(u(1)))/2;
    else
        ret = trapezoidalIntegral(dx,absSq(u));
    end
end