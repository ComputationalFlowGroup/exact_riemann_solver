function [fofpprime] = fofpprime(pguess,g,WL,WR)
%FOFPPRIME SummARy of this function goes here
%   Detailed expLanation goes here

rhoL = WL(1); uL = WL(2); pL = WL(3);

rhoR = WR(1); uR = WR(2); pR = WR(3);


AL = 2/((g+1)*rhoL); % shock pARameter
BL = ((g-1)/(g+1))*pL; % shock pARameter
cL = sqrt(g*pL/rhoL); % speed of sound

AR = 2/((g+1)*rhoR); % shock pARameter
BR = ((g-1)/(g+1))*pR; % shock pARameter
cR = sqrt(g*pR/rhoR); % speed of sound


if pguess > pL % shock state, but it is now fL'(p) not fL(p)
    fLp = ((AL/(BL+pguess))^(1/2))*(1-0.5*(pguess-pL)/(pguess+BL)); % left shock
else % rarefaction case
    fLp = (1/(rhoL*cL))*(pguess/pL)^(-(g+1)/(2*g)); % left rarefaction
end

if pguess > pR % shock state, but it is now fR'(p) not fR(p)
    fRp = ((AR/(BR+pguess))^(1/2))*(1-0.5*(pguess-pR)/(pguess+BR)); % left shock
else % rarefaction case
    fRp = (1/(rhoR*cR))*(pguess/pR)^(-(g+1)/(2*g)); % right rarefaction
end

fofpprime = fLp + fRp; 

end