function spDeriv = finiteDifferenceSpMatrix(derivPoints,n,M)
% finiteDifferenceSpMatrix  - Returns a square, sparse matrix which, when
%                             multiplied with a vector containing function
%                             values, returns a derivative approximation
% Syntax: spDeriv = finiteDifferenceSpMatrix(derivPoints,n,M)
%
% Input:
% derivPoints   - How many points, including the point we are approximating
%                 the derivative in, are used to approximate the derivative.
% derivNo       - The order of the derivative.
% M             - The size of the matrix, MxM.
%
% Output:
% spDeriv   - A sparse MxM matrix with proper finite difference
%             coefficients in order to approximate the n'th derivative.
%
% Non-standard dependencies: None.
% See also: initModelInfo.m, finiteDifferenceSpMatrix.m.
%
% Note: 
% derivPoints odd required
% n < derivPoints required
% derivPoints <= M required

if mod(derivPoints,2)==0
    error('Centered derivative approximation required. derivPoints odd required.')
elseif n >= derivPoints
    error('Derivative approximations requires more points than the order. n < derivPoints required. ')
elseif derivPoints > M
    error('Number of points used to approximate derivative cannot exceed size of matrix. derivPoints <= M required. ')
end

% Calculate the vector
vector = finiteDifferenceVector(derivPoints,n);

% Since it's periodic, we have that some of the steps will wrap around
noWrap = (derivPoints-1)/2;

diagMat = zeros(M,derivPoints);
for i = 1:derivPoints
    diagMat(:,i) = vector(i)*ones(M,1);
end

% Construct the corner elements
cornerCoordVec1 = zeros(1,noWrap*(noWrap+1));
valVec = zeros(1,noWrap*(noWrap+1));

currFilledCoord = 0;
currFilledVal = 0;
for i = 1:noWrap
    % Create the coordinates
    cornerCoordVec1((currFilledCoord+1):(currFilledCoord+1+(noWrap-i))) = 1:(noWrap-i+1);
    cornerCoordVec2((currFilledCoord+1):(currFilledCoord+1+(noWrap-i))) = (M-noWrap+i):M;
    currFilledCoord = currFilledCoord + noWrap-i+1;
    cornerCoordVec1((currFilledCoord+1):(currFilledCoord+1+(noWrap-i))) = (M-noWrap+i):M;
    cornerCoordVec2((currFilledCoord+1):(currFilledCoord+1+(noWrap-i))) = 1:(noWrap-i+1);
    currFilledCoord = currFilledCoord + noWrap-i+1;
    % Fill the values
    valVec((currFilledVal+1):(currFilledVal+2*(noWrap-i+1))) = vector(i)*ones(1,2*(noWrap-i+1));
    currFilledVal = currFilledVal + 2*(noWrap-i+1);
end

% Now add the diagonals with the corner elements
spDeriv = spdiags(diagMat,-noWrap:noWrap,M,M) +...
    sparse(cornerCoordVec1,cornerCoordVec2,valVec);
end