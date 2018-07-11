function weight = finiteDifferenceVector(derivPoints,n)
% finiteDifferenceVector    - Calculates how to weight the points in order 
%                             to approximate the n'th derivative.
% Syntax: weight = finiteDifferenceVector(derivPoints,n)
%
% Input:
% derivPoints   - How many points, including the point we are approximating
%                 the derivative in, are used to approximate the derivative.
% n             - The order of the derivative.
%
% Output:
% weight    - A vector of length derivPoints containing the weights.
%
% Non-standard dependencies: None.
% See also: initModelInfo.m, finiteDifferenceSpMatrix.m.
%
% Made with the help of http://web.media.mit.edu/~crtaylor/calculator.html

bound = (derivPoints-1)/2;

row = -bound:bound;

mat = zeros(derivPoints);

for i = 1:derivPoints
    mat(i,:) = row.^(i-1);
end

vec = zeros(derivPoints,1);
vec(n+1) = factorial(n);

weight = mat\vec; % /dx^2

end