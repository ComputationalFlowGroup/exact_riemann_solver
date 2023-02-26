function [pstar] = euler_riemann_root_finder(pguess,WL,WR,g,tolerance)
%EULER_RIEMANN_ROOT_FINDER Summary of this function goes here
%   Solving for the pstar and ustar of the Riemann problem to the Euler
%   equations


error = 1; % error to start the while loop

% what we want to find is 
% f(pstar) = 0, this is the equation we want to solve
% f(pguess) = c, we will use this to define an error
% error = abs(f(pguess) - f(pnewguess)), we are trying to drive this error
% to zero
% pnewguess can be calculated by using a Taylor series expansion of f(p+dp)

% the Taylor series expansion of f(p+dp)

% f(p+dp) = f(p) + f'(p)*dp + O(dp^2)

% f(pguess+dp) = f(pguess) + f'(pguess)*(pnewguess-pguess), if we arrived at our solution
% f(pguess+dp) = 0, this is such that dp = 0 and p = pstar

% 0 = f(pguess) + f'(pguess)*(pnewguess-pguess), rearranging, we can find
% an equation for pnewguess

% pnewguess = pguess - f(pguess)/f'(pguess), we will use this to then find
% pstar

while error > tolerance % this is a while loop to keep iterating to the solution
    pnewguess = pguess - fofp(pguess,g,WL,WR)/fofpprime(pguess,g,WL,WR);
    %error = abs(fofp(pguess) - fofp(pnewguess));
    % calculating the error
    error = abs((pnewguess-pguess))/(0.5*abs((pguess+pnewguess)));
    % updating the pguess, to continue iterating
    pguess = pnewguess;
end
rhoL = WL(1); uL = WL(2); pL = WL(3);
rhoR = WR(1); uR = WR(2); pR = WR(3);
if pguess < 0 && uL < 0 && uR > 0
    aL = sqrt(g*pL/rhoL); % speed of sound
    aR = sqrt(g*pR/rhoR); % speed of sound    
    num = aL + aR - 0.5*(g-1)*(uR-uL);
    dem = aL/(pL)^((g-1)/(2*g)) + aR/(pR)^((g-1)/(2*g));
    pguess = (num/dem)^(2*g/(g-1));
end

% out of this while loop finished iterating
pstar = pguess;

end