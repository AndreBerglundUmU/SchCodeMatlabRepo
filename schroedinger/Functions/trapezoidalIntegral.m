function [integral] = trapezoidalIntegral(dx,f)
% trapezoidalIntegral   - The trapezoidal integral approximation of f
% Syntax: [integral] = trapezoidalIntegral(dx,f)
%
% Input:
% dx    - Space step size
% f     - A vector of f(x), where x is a vector such that x(i)-x(i-1)=dx
%
% Output:
% integral  - A scalar of the trapezoidal integral approximation
%
% Non-standard dependencies: None.
% See also: Any accompanying script for example usage.
integral = dx/2*sum((f(1:(end-1)) + f(2:end)));
end

