function [f] = f_ps(p,g,WL,WR)
%RIEMANN_SOLVER Summary of this function goes here
%   Detailed explanation goes here
rl = WL(1); ul = WL(2); pl = WL(3);
rr = WR(1); ur = WR(2); pr = WR(3);
al = 2/((g+1)*rl); bl = pl*(g-1)/(g+1);
cl = sqrt(g*pl/rl);
ar = 2/((g+1)*rr); br = pr*(g-1)/(g+1);
cr = sqrt(g*pr/rr);
if p > pl
    fl = (p-pl).*(al./(p+bl)).^0.5; % left shock
else
    fl = (2*cl/(g-1))*((p./pl).^((g-1)/(2*g))-1); % left rarefaction
end
if p > pr
    fr = (p-pr).*(ar./(p+br)).^0.5; % left shock
else
    fr = (2*cr/(g-1))*((p./pr).^((g-1)/(2*g))-1); % left rarefaction
end
f = fl + fr + ur - ul;
end