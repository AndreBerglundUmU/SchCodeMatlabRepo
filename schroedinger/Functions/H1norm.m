function ret = H1norm(u,dx,k,per)
% H1norm - A trapezoidal integral approximation of the H1 norm.
% Syntax: ret = H1norm(u,dx,k,per)
%
% Input:
% u     - A vector of length M containing the function u in Fourier space.
% dx    - The space step size.
% k     - A vector of length M containing the Fouier modes.
% per   - A boolean value declaring whether u periodic or not.
%
% Output:
% ret - int abs(u)^2 + abs(d/dx u)^2 dx
%
% Non-standard dependencies: L2norm.m.
% See also: L2norm.m
    ret = L2norm(ifft(u),dx,per) + L2norm(ifft(1i*k.*u),dx,per);
end