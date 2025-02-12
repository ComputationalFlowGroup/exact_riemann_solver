function [fL] = fLofp(pguess,g,WL)
%FLOFP Summary of this function goes here
%   Detailed explanation goes here

% we need to account for both the shock and the rarefaction
% let's first calculate a few things with WL

rhoL = WL(1); uL = WL(2); pL = WL(3); % primitive variables
A = 2/((g+1)*rhoL); % shock parameter
B = ((g-1)/(g+1))*pL; % shock parameter
cl = sqrt(g*pL/rhoL); % speed of sound

if pguess > pL % left shock 
    fL = (pguess - pL) * sqrt(A/(pguess+B));
else % left rarefaction
    fL = 2*cl/(g-1) * ((pguess/pL).^((g-1)/(2*g)) - 1);
end

end