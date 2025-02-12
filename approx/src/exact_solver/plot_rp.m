function [rsol] = plot_rp(ps, g, WL, WR, LX0, LX1, dloc, N, T)
%PLOT_RP Summary of this function goes here
%   Detailed explanation goes here
mu = (g-1)/(g+1);
rl = WL(1); ul = WL(2); pl = WL(3);
rr = WR(1); ur = WR(2); pr = WR(3);
Al = 2/((g+1)*rl); Bl = pl*mu;
cl = sqrt(g*pl/rl);
Ar = 2/((g+1)*rr); Br = pr*mu;
cr = sqrt(g*pr/rr);

dx = (LX1-LX0)/N;
x = ((LX0+dx/2):dx:(LX1-dx/2));
x = x-dloc;

rsol = zeros(N,1); usol = zeros(N,1); psol = zeros(N,1);

for i = 1:N
    %left state
    if ps > pl      % left shock
        rsl = rl*((ps/pl + mu)/(mu*(ps/pl) + 1));
        usl = ul - (ps-pl)*(Al/(ps+Bl))^0.5;
        Sl = ul - cl*(((g+1)/(2*g))*(ps/pl) + (g-1)/(2*g))^0.5;
        if x(i) < Sl*T
            rsol(i) = rl;
            usol(i) = ul;
            psol(i) = pl;
        elseif (x(i) > Sl*T && x(i) < usl*T)
            rsol(i) = rsl;
            usol(i) = usl;
            psol(i) = ps;
        end        
    else            % left rarefaction
        rsl = rl*(ps/pl)^(1/g);
        usl = ul - ((2*cl)/(g-1))*((ps/pl)^((g-1)/(2*g)) - 1);        
        csl = cl*(ps/pl)^((g-1)/(2*g));        
        Shl = ul - cl;        
        Stl = usl - csl;
        if (x(i) < Shl*T)
            rsol(i) = rl;
            usol(i) = ul;
            psol(i) = pl;
        elseif (x(i) > Shl*T && x(i) < Stl*T)
            rsol(i) = rl*(2/(g+1) + (mu/cl)*(ul-x(i)/T))^(2/(g-1));    
            usol(i) = (2/(g+1))*(cl+ul*(g-1)/2 + x(i)/T);
            psol(i) = pl*(2/(g+1) + (mu/cl)*(ul-x(i)/T))^((2*g)/(g-1));
        elseif (x(i) > Stl*T && x(i) < usl*T)
            rsol(i) = rsl;
            usol(i) = usl;
            psol(i) = ps;
        end           
    end
    
    %right state
    if ps > pr      % right shock
        rsr = rr*((ps/pr + mu)/(mu*(ps/pr) + 1));
        usr = ur + (ps-pr)*(Ar/(ps+Br))^0.5;
        Sr = ur + cr*(((g+1)/(2*g))*(ps/pr) + (g-1)/(2*g))^0.5;
        if (x(i) >= Sr*T)
            rsol(i) = rr;
            usol(i) = ur;
            psol(i) = pr;
        elseif (x(i) <= Sr*T && x(i) >= usr*T)
            rsol(i) = rsr;
            usol(i) = usr;
            psol(i) = ps;
        end        
    else            % right rarefaction
        rsr = rr*(ps/pr)^(1/g);
        usr = ur + ((2*cr)/(g-1))*((ps/pr)^((g-1)/(2*g)) - 1);        
        csr = cr*(ps/pr)^((g-1)/(2*g));        
        Shr = ur + cr;        
        Str = usr + csr;
        if (x(i) > Shr*T)
            rsol(i) = rr;
            usol(i) = ur;
            psol(i) = pr;
        elseif (x(i) < Shr*T && x(i) > Str*T)
            rsol(i) = rr*(2/(g+1) - (mu/cr)*(ur-x(i)/T))^(2/(g-1));    
            usol(i) = (2/(g+1))*(-cr+ur*(g-1)/2 + x(i)/T);
            psol(i) = pr*(2/(g+1)-(mu/cr)*(ur-x(i)/T))^((2*g)/(g-1));
        elseif (x(i) < Str*T && x(i) > usr*T)
            rsol(i) = rsr;
            usol(i) = usr;
            psol(i) = ps;
        end           
    end
   
end

x = x + dloc;
esol = psol./((g-1).*rsol);
subplot(2,2,1)
hold on;
plot(x,rsol,'k')
ylabel('$\rho$','Interpreter','Latex','FontSize',14);
set(gca,'FontName','Times','FontSize',14);
subplot(2,2,2)
hold on;
plot(x,usol,'k')
ylabel('u','Interpreter','Latex','FontSize',14);
set(gca,'FontName','Times','FontSize',14);
subplot(2,2,3)
hold on;
plot(x,psol,'k')
ylabel('p','Interpreter','Latex','FontSize',14);
set(gca,'FontName','Times','FontSize',14);
subplot(2,2,4)
hold on;
plot(x,esol,'k')
ylabel('internal e','Interpreter','Latex','FontSize',14);
set(gca,'FontName','Times','FontSize',14);
end