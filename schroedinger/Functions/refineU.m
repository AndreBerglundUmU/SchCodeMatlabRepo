function ret = refineU(u)
    ret = [u ; (u + [u(2:end) u(1)])/2];
    ret = ret(:)';
end