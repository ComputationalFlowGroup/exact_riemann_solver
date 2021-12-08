function [f] = fp_ps(p,g,WL,WR)
%RIEMANN_SOLVER Summary of this function goes here
%   Detailed explanation goes here
rl = WL(1); ul = WL(2); pl = WL(3);
rr = WR(1); ur = WR(2); pr = WR(3);
al = 2/((g+1)*rl); bl = pl*(g-1)/(g+1);
cl = sqrt(g*pl/rl);
ar = 2/((g+1)*rr); br = pr*(g-1)/(g+1);
cr = sqrt(g*pr/rr);

if p > pl
    fl = ((al/(bl+p))^(1/2))*(1-(p-pl)/(2*(bl+p))); % left shock
else
    fl = (1/(rl*cl))*(p/pl)^(-(g-1)/(2*g)); % left rarefaction
end

if p > pr
    fr = ((ar/(br+p))^(1/2))*(1-(p-pr)/(2*(br+p))); % right shock
else
    fr = (1/(rr*cr))*(p/pr)^(-(g-1)/(2*g)); % right rarefaction
end

f = fl + fr;

end