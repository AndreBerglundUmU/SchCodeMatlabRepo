function weight = finiteDifferenceVector(derivPoints,derivNo)
% Odd number of points assumed, gives even weighting
% DerivNo < derivPoints required

bound = (derivPoints-1)/2;

row = -bound:bound;

mat = zeros(derivPoints);

for i = 1:derivPoints
    mat(i,:) = row.^(i-1);
end

vec = zeros(derivPoints,1);
vec(derivNo+1) = factorial(derivNo);

weight = mat\vec; % /dx^2

end