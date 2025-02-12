function [] = rp_plot(wL, wsL, wR, wsR, LX0, LX1, dloc, N, T, filename, solname,sty)
%PLOT_RP Summary of this function goes here
%   Detailed explanation goes here

% left state variables
gL = wL(1); cL = wL(2); bL = wL(3); rL = wL(4); uL = wL(5); vL = wL(6);
pL = wL(7); piL= wL(8); s11L = wL(9); s12L = wL(10); s22L = wL(11);
csL = sqrt( cL*cL + (4/3)*bL/(rL*rL) );
shL = uL - csL;

% left star state variables
rsL = wsL(4); usL = wsL(5); psL = wsL(7); s11sL = wsL(9); s22sL = wsL(11);
cssL = sqrt( cL*cL*((psL+piL)/(pL+piL))^((gL-1)/gL) + (4/3)*bL/(rsL*rsL) );
stL = usL - cssL;
mL = sqrt( rL*rL*cL*cL*(((gL+1)/(2*gL))*((psL+piL)/(pL+piL))+(gL-1)/(2*gL))+(4/3)*bL );
sL = usL - mL/rsL;

% right state variables
gR = wR(1); cR = wR(2); bR = wR(3); rR = wR(4); uR = wR(5); vR = wR(6);
pR = wR(7); piR= wR(8); s11R = wR(9); s12R = wR(10); s22R = wR(11);
csR = sqrt( cR*cR + (4/3)*bR/(rR*rR) );
stR = uR + csR;

% right star state variables
rsR = wsR(4); usR = wsR(5); psR = wsR(7); s11sR = wsR(9); s22sR = wsR(11);
cssR = sqrt( cR*cR*((psR+piR)/(pR+piR))^((gR-1)/gR) + (4/3)*bR/(rsR*rsR) );
shR = usR + cssR;
mR = sqrt( rR*rR*cR*cR*(((gR+1)/(2*gR))*((psR+piR)/(pR+piR))+(gR-1)/(2*gR))+(4/3)*bR );
sR = usR + mR/rsR;

% shear star state
% ssL = 0.5*(uL+usL) - sqrt(bL/(rL*rL));
ssL = usL - sqrt(bL/(rL*rL));
% ssL = uL - mL/rL;
% ssR = 0.5*(uR+usR) + sqrt(bR/(rR*rR));
ssR = usR + sqrt(bR/(rR*rR));
% ssR = uR + mR/rR;
s12s = (s12R + s12L + sqrt(bR)*(vR-vL)) / 2;
vs = (sqrt(bR)*vL+sqrt(bR)*vR+s12R-s12L) / (2*sqrt(bR));

% grid calculation
dx = (LX1-LX0)/N; x = ((LX0+dx/2):dx:(LX1-dx/2)); x = x-dloc;

% solution variables
rsol = zeros(N,1); usol = zeros(N,1); psol = zeros(N,1); 
s11sol = zeros(N,1); vsol = zeros(N,1); s12sol = zeros(N,1);
s22sol = zeros(N,1);

for i = 1:N
    %left state
    if psL > pL      % left shock
        if x(i) <= sL*T
            rsol(i) = rL;
            usol(i) = uL;
            psol(i) = pL;
            s11sol(i) = s11L;
            s22sol(i) = s22L;
        elseif (x(i) > sL*T && x(i) < usL*T)
            rsol(i) = rsL;
            usol(i) = usL;
            psol(i) = psL;
            s11sol(i) = s11sL;
            s22sol(i) = s22sL;
        end        
    else            % left rarefaction
        if (x(i) < shL*T)
            rsol(i) = rL;
            usol(i) = uL;
            psol(i) = pL;
            s11sol(i) = s11L;
            s22sol(i) = s22L;
        elseif (x(i) > shL*T && x(i) < stL*T)
            [r,u,p,s11,s22] = ps_rare_calc(x(i)/T,wL,-1);
            rsol(i) = r;    
            usol(i) = u;
            psol(i) = p;
            s11sol(i) = s11;
            s22sol(i) = s22;
        elseif (x(i) > stL*T && x(i) < usL*T)
            rsol(i) = rsL;
            usol(i) = usL;
            psol(i) = psL;
            s11sol(i) = s11sL;
            s22sol(i) = s22sL;
        end           
    end
    
    %right state
    if psR > pR      % right shock
        if (x(i) >= sR*T)
            rsol(i) = rR;
            usol(i) = uR;
            psol(i) = pR;
            s11sol(i) = s11R;
            s22sol(i) = s22R;
        elseif (x(i) < sR*T && x(i) >= usR*T)
            rsol(i) = rsR;
            usol(i) = usR;
            psol(i) = psR;
            s11sol(i) = s11sR;
            s22sol(i) = s22sR;
        end
    else            % right rarefaction
        if (x(i) > shR*T)
            rsol(i) = rR;
            usol(i) = uR;
            psol(i) = pR;
            s11sol(i) = s11R;
            s22sol(i) = s22R;
        elseif (x(i) < shR*T && x(i) > stR*T)
            [r,u,p,s11,s22] = ps_rare_calc(x(i)/T,wR,1);
            rsol(i) = r;    
            usol(i) = u;
            psol(i) = p;
            s11sol(i) = s11;
            s22sol(i) = s22;
        elseif (x(i) < stR*T && x(i) > usR*T)
            rsol(i) = rsR;
            usol(i) = usR;
            psol(i) = psR;
            s11sol(i) = s11sR;
            s22sol(i) = s22sR;
        end           
    end
    
    % shear stress
    if x(i) < ssL*T
        vsol(i) = vL;
        s12sol(i) = s12L;
    elseif (x(i) < ssR*T && x(i) > ssL*T) % middle state
        vsol(i) = vs;
        s12sol(i) = s12s;
    elseif (x(i) > ssR*T)
        vsol(i) = vR;
        s12sol(i) = s12R;
    end
end

x = x + dloc;
nsol = dlmread(solname, ' ', 1, 0);
figure(1)
subplot(2,3,1)
hold on;
plot(x,rsol,'-k');
plot(nsol(:,1),nsol(:,2),sty);
ylabel('$\rho$','Interpreter','Latex','FontSize',14);
subplot(2,3,2)
hold on;
plot(x,usol,'-k');
plot(nsol(:,1),nsol(:,3),sty);
ylabel('$u$','Interpreter','Latex','FontSize',14);
subplot(2,3,3)
hold on;
plot(x,psol,'-k');
plot(nsol(:,1),nsol(:,5),sty);
ylabel('$p$','Interpreter','Latex','FontSize',14);
subplot(2,3,4)
hold on;
plot(x,(psol - s11sol)/1e6,'-k');
plot(nsol(:,1),(nsol(:,5)-nsol(:,7))/1e6,sty);
% plot(x,s11sol,'-k');
% plot(nsol(:,1),nsol(:,7),sty);
ylabel('$\sigma_{11}$','Interpreter','Latex','FontSize',14);
subplot(2,3,5)
hold on;
plot(x,s12sol,'-k');
plot(nsol(:,1),nsol(:,9),sty);
ylabel('$\tau_{12}$','Interpreter','Latex','FontSize',14);
subplot(2,3,6)
hold on;
plot(x,vsol,'-k');
plot(nsol(:,1),nsol(:,4),sty);
ylabel('$-\sigma_{22}$','Interpreter','Latex','FontSize',14);
% set(gca,'FontName','Times','FontSize',14);
A = [x',rsol,usol,vsol,psol,s11sol,s22sol,s12sol];
fileID = fopen(filename,'w');
fprintf(fileID,'%10s\n','VARIABLES="X","RHO","U","V","P","S11","S22","S12"');
fprintf(fileID,'%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n',A');
fclose(fileID);
end

% rarefaction solver
function [rs,us,ps,s11s,s22s] = ps_rare_calc(eig,w0,sign)
    g0 = w0(1);
    c0 = w0(2);
    b0 = w0(3);
    r0 = w0(4);
    u0 = w0(5);
    p0 = w0(7);
    pi0= w0(8); 
    s110 = w0(9);
    s220 = w0(11);
    K0 = ((p0+pi0)^(1/g0))/(r0*g0);
    mu0 = (g0+1)/g0;
    fps = @(p) sqrt( K0*p.^(-mu0) + (4*b0/3)*K0*K0*p.^(-2*mu0) );
    cps = @(p) sqrt( c0*c0*((p+pi0)/(p0+pi0)).^((g0-1)/g0) +...
        ((4*b0)/(3.*r0.*r0))*((p+pi0)/(p0+pi0)).^(-2/g0) );
    ifps = @(ps) u0 - (eig - sign*cps(ps)) + sign*integral(fps,p0+pi0,ps+pi0) ;
    ps = fzero(ifps,p0);
    us = eig - sign*cps(ps);
    xi = ((ps+pi0)/(p0+pi0))^(-1/g0);
    s11s = s110 - (4*b0/(3*r0))*(1-xi);
    xi = ((ps+pi0)/(p0+pi0))^(1/g0);
    rs = r0*xi;
    s22s = s220 - (2*b0/(3*r0))*(1-xi);
end