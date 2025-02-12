function [ flx ] = roe_flux_e(U,l,r,g,dx,dt,flim)
%ROE_FLUX Summary of this function goes here
%   Detailed explanation goes here

pv = @(U,j) (g-1)*(U(j,3)-0.5*(U(j,2)*U(j,2)/U(j,1))); 

flx_vec = @(U,j) [U(j,2),(U(j,2)*U(j,2)/U(j,1))+pv(U,j),...
    (U(j,2)/U(j,1))*(U(j,3)+pv(U,j))];

rl = U(l,1);
rr = U(r,1);
rlroe = sqrt(rl);
rrroe = sqrt(rr);
rul = U(l,2);
rur = U(r,2);
ul = U(l,2)/rl;
al = sqrt(g*pv(U,l)/rl);
laml = ul-al;
ur = U(r,2)/rr;
ar = sqrt(g*pv(U,r)/rr);
lamr = ur+ar;
hl = (U(l,3)+pv(U,l))/rl;
hr = (U(r,3)+pv(U,r))/rr;

uroe = (rlroe*ul+rrroe*ur)/(rlroe+rrroe);
hroe = (rlroe*hl+rrroe*hr)/(rlroe+rrroe);
aroe = sqrt((g-1)*(hroe-0.5*uroe*uroe));

al2 = ((g-1)/(aroe*aroe))*((rr-rl)*(hroe-uroe*uroe)+uroe*(rur-rul)-(U(r,3)-U(l,3)));
al1 = (1/(2*aroe))*((rr-rl)*(uroe+aroe) - (rur-rul) - aroe*al2);
al3 = (rr-rl)-(al1+al2);

alpha = [al1, al2, al3];  
lam = [uroe-aroe, uroe, uroe+aroe];
K = [1, uroe-aroe, hroe-uroe*aroe;...
      1, uroe, 0.5*uroe*uroe;...
      1, uroe+aroe, hroe+uroe*aroe];  

ustar = (rl*ul+al1*(uroe-aroe))/(rl+al1);
rlstar = rl + al1;
rrstar = rr - al3;
pstar = (g-1)*(U(l,3)+al1*(hroe-uroe*aroe)-0.5*(rlstar*ustar*ustar));
alstar = sqrt(g*pstar/rlstar);
arstar = sqrt(g*pstar/rrstar);
lamls = ustar-alstar;
lamrs = ustar+arstar;
lam1b = laml*((lamls-lam(1))/(lamls-laml));
lam3b = lamr*((lam(3)-lamrs)/(lamr-lamrs));

fl = flx_vec(U,l);
diff = zeros(1,3);
nu = dt/dx*lam;

pj = l-sign(nu);
phi = zeros(1,3);

for j = 1:3
    phi(j) = flx_lim(U, pj, j, g, alpha, flim);
end

if ((laml < 0 && lamls > 0)) 
    flx = fl+lam1b*al1*K(1,:) + 0.5*lam1b*(sign(nu(1))-nu(1))*al1*K(1,:)*phi(1);
elseif (lamrs < 0 && lamr > 0)
    flx = fl-lam3b*al3*K(3,:) + 0.5*lam3b*(sign(nu(3))-nu(3))*al3*K(3,:)*phi(3);
else
    for ind = 1:3
        diff = diff + alpha(ind).*min([lam(ind),0]).*K(ind,:)...
            + 0.5*lam(ind)*(sign(nu(ind))-nu(ind))*alpha(ind)*K(ind,:)*phi(ind);
    end
    flx = fl+diff;
end
   
end