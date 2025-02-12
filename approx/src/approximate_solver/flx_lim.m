function [ phi ] = flx_lim( U, pj, j, g, alpha, flim)
%FLX_LIM Summary of this function goes here
%   Detailed explanation goes here

l = pj(j); 
r = l + 1;

pv = @(U,j) (g-1)*(U(j,3)-0.5*(U(j,2)*U(j,2)/U(j,1))); 

rl = U(l,1);
rr = U(r,1);
rlroe = sqrt(rl);
rrroe = sqrt(rr);
rul = U(l,2);
rur = U(r,2);
ul = U(l,2)/rl;
ur = U(r,2)/rr;
hl = (U(l,3)+pv(U,l))/rl;
hr = (U(r,3)+pv(U,r))/rr;

uroe = (rlroe*ul+rrroe*ur)/(rlroe+rrroe);
hroe = (rlroe*hl+rrroe*hr)/(rlroe+rrroe);
aroe = sqrt((g-1)*(hroe-0.5*uroe*uroe));

al2 = ((g-1)/(aroe*aroe))*((rr-rl)*(hroe-uroe*uroe)+uroe*(rur-rul)-(U(r,3)-U(l,3)));
al1 = (1/(2*aroe))*((rr-rl)*(uroe+aroe) - (rur-rul) - aroe*al2);
al3 = (rr-rl)-(al1+al2);
alp = [al1, al2, al3];

theta = alp(j)/alpha(j);

switch (flim)
    case 0 % superbee
        phi = max([0,min([1,2*theta]),min([theta,2])]);
    case 1 % minmod 
        phi = max([0,min([1,theta])]);
end

end