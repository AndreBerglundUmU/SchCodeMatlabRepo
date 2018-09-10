function ret = refineU(u)
% refineU   - Refine u using linear interpolation. Assuming u periodic.
% Syntax: ret = refineU(u)
%
% Input:
% u     - A vector of size (1,M) containing the function u
%
% Output:
% ret   - A vector of size (1,2*M) containining the refined u
%
% Non-standard dependencies: None.
% See also: ~.
    if size(u,1)~= 1
        error('Error: u must be a vector of size (1,M).')
    end
    ret = [u ; (u + [u(2:end) u(1)])/2];
    ret = ret(:)';
end


% 
% function ret = periodicRefineVector(x)
%     ret = zeros(1,2*length(x));
%     ret(1:2:(end-1)) = x;
%     ret(2:2:end) = ([x(2:end),x(1)]-x(1:end))/2;
% end