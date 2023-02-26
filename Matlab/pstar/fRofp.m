function [fR] = fRofp(pguess,g,WR)
%FROFP Summary of this function goes here
%   Detailed explanation goes here

% we need to account for both the shock and the rarefaction
% let's first calculate a few things with WR

rhoR = WR(1); uR = WR(2); pR = WR(3); % primitive variables
A = 2/((g+1)*rhoR); % shock parameter
B = ((g-1)/(g+1))*pR; % shock parameter
cl = sqrt(g*pR/rhoR); % speed of sound

if pguess > pR % left shock 
    fR = (pguess - pR) * sqrt(A/(pguess+B));
else % left rarefaction
    fR = 2*cl/(g-1) * ((pguess/pR).^((g-1)/(2*g)) - 1);
end

end

