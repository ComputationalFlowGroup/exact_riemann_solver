function [rsol] = plot_roe_eue(g,WL,WR,LX0,LX1,NX,T,cfl,flim)
%LXF Summary of this function goes here
%   Detailed explanation goes here

uvl = @(W) [W(1),W(1)*W(2),W(3)/(g-1)+0.5*W(1)*W(2)*W(2)];

dx = (LX1-LX0)/NX;

x = (LX0+dx/2):dx:(LX1-dx/2);
usl = zeros(NX,3);

for ind = 1:NX
    if ind < NX/2+1
        usl(ind,:) = uvl(WL);
    else
        usl(ind,:) = uvl(WR);
    end
end
utemp = usl;

t = 0;
a_calc = @(usl) usl(:,2)./usl(:,1) + sqrt(g*((g-1)*(usl(:,3)-0.5*usl(:,2).*usl(:,2)./usl(:,1)))./usl(:,1));
dt_calc = @(usl) cfl*dx/max(abs(a_calc(usl)));
dt = dt_calc(usl);

while t < T
    flxl = roe_flux_e(usl,2,3,g,dx,dt,flim);
    for ind = 3:NX-2
        flxr = roe_flux_e(usl,ind,ind+1,g,dx,dt,flim);
        utemp(ind,:) = usl(ind,:) - (dt/dx)*(flxr-flxl);
        flxl = flxr;
    end
    usl = utemp;
    t = t + dt;
    dt = dt_calc(usl);
end

rsol = usl(:,1);
usol = usl(:,2)./usl(:,1);
psol = (g-1)*(usl(:,3)-0.5*usl(:,2).*usl(:,2)./usl(:,1));
esol = psol./((g-1).*rsol);

subplot(2,2,1)
hold on;
plot(x,rsol,'.k')
ylabel('$\rho$','Interpreter','Latex','FontSize',14);
set(gca,'FontName','Times','FontSize',14);
subplot(2,2,2)
hold on;
plot(x,usol,'.k')
ylabel('u','Interpreter','Latex','FontSize',14);
set(gca,'FontName','Times','FontSize',14);
subplot(2,2,3)
hold on;
plot(x,psol,'.k')
ylabel('p','Interpreter','Latex','FontSize',14);
set(gca,'FontName','Times','FontSize',14);
subplot(2,2,4)
hold on;
plot(x,esol,'.k')
ylabel('internal e','Interpreter','Latex','FontSize',14);
set(gca,'FontName','Times','FontSize',14);

end