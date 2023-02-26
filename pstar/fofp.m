function [f] = fofp(pguess,g,WL,WR)
%FOFP Summary of this function goes here
%   Detailed explanation goes here

uL = WL(2);
uR = WR(2);

fL = fLofp(pguess,g,WL);
fR = fRofp(pguess,g,WR);
f =  fL + fR  + uR - uL;

end