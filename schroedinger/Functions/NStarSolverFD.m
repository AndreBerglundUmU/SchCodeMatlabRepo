function [ NStar ] = NStarSolverFD(currU,a,h,sigma)
crit = true;
K = length(currU);
NStar = currU;
i = 1;
while crit && i<120
    oldNStar = NStar;
    tempNStar = a+h/2*NStar;
    NStar = 1i*nonLin(tempNStar,sigma);
    crit = norm((oldNStar-NStar)./K,2) > eps;
    i = i+1;
end
end

