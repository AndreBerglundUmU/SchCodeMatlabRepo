function [integral] = trapezoidalIntegral(dx,f)
% This function will approximate the integral int(f(x)) given the input:
%   dx  - Space step size
%   f   - A vector of f(x), where x is a vector such that x(i)-x(i-1)=dx
integral = dx/2*sum((f(1:(end-1)) + f(2:end)));
end

