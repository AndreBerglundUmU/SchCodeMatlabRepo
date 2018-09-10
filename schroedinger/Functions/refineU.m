function ret = refineU(u,periodic)
% refineU   - Refine u using linear interpolation. Assuming u periodic.
% Syntax: ret = refineU(u)
%
% Input:
% u         - A vector of size (1,M) containing the function u
% periodic  - A boolean value which decides if u is interpreted as periodic
%
% Output:
% ret   - A vector containining the refined u, of size (1,2*M) if periodic
%         or of size (1,2*M-1) if not
%
% Non-standard dependencies: None.
% See also: ~.
    if size(u,1)~= 1
        error('Error: u must be a vector of size (1,M).')
    end
    
    if periodic
        ret = [u ; (u + [u(2:end) u(1)])/2];
        ret = ret(:)';
    else
        M = length(u);
        ret = zeros(1,2*M-1);
        ret(1:2:end) = u;
        ret(2:2:end-1) = (u(1,end-1) + u(2:end))/2;
    end
end


% 
% function ret = periodicRefineVector(x)
%     ret = zeros(1,2*length(x));
%     ret(1:2:(end-1)) = x;
%     ret(2:2:end) = ([x(2:end),x(1)]-x(1:end))/2;
% end