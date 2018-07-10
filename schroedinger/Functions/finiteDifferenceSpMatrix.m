function Delta = finiteDifferenceSpMatrix(derivPoints,derivNo,noPoints)
% This will create a periodic finite difference matrix
% derivPoints odd required
% derivNo < derivPoints required
% derivPoints <= noPoints required

% Calculate the vector
vector = finiteDifferenceVector(derivPoints,derivNo);

% Since it's periodic, we have that some of the steps will wrap around
noWrap = (derivPoints-1)/2;

diagMat = zeros(noPoints,derivPoints);
for i = 1:derivPoints
    diagMat(:,i) = vector(i)*ones(noPoints,1);
end

% Construct the corner elements
cornerCoordVec1 = zeros(1,noWrap*(noWrap+1));
valVec = zeros(1,noWrap*(noWrap+1));

currFilledCoord = 0;
currFilledVal = 0;
for i = 1:noWrap
    % Create the coordinates
    cornerCoordVec1((currFilledCoord+1):(currFilledCoord+1+(noWrap-i))) = 1:(noWrap-i+1);
    cornerCoordVec2((currFilledCoord+1):(currFilledCoord+1+(noWrap-i))) = (noPoints-noWrap+i):noPoints;
    currFilledCoord = currFilledCoord + noWrap-i+1;
    cornerCoordVec1((currFilledCoord+1):(currFilledCoord+1+(noWrap-i))) = (noPoints-noWrap+i):noPoints;
    cornerCoordVec2((currFilledCoord+1):(currFilledCoord+1+(noWrap-i))) = 1:(noWrap-i+1);
    currFilledCoord = currFilledCoord + noWrap-i+1;
    % Create the values
    valVec((currFilledVal+1):(currFilledVal+2*(noWrap-i+1))) = vector(i)*ones(1,2*(noWrap-i+1));
    currFilledVal = currFilledVal + 2*(noWrap-i+1);
end

% Now add the diagonals with the corner elements
Delta = spdiags(diagMat,-noWrap:noWrap,noPoints,noPoints) +...
    sparse(cornerCoordVec1,cornerCoordVec2,valVec);
end