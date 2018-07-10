function ret = H1norm(u,dx,k,per)
    ret = L2norm(ifft(u),dx,per) + L2norm(ifft(1i*k.*u),dx,per);
end