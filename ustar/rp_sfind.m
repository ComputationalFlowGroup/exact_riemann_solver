function [ wL, wR ] = rp_sfind( wL, wR )
%PSTAR_FIND Summary of this function goes here
%   Detailed explanation goes here

% Left state
cL = wL(2); bL = wL(3); rL = wL(4); uL = wL(5); pL = wL(7); s11L = wL(9);
zL = rL*cL*cL + (4/3)*bL/(rL*rL);

% Right state
cR = wR(2); bR = wR(3); rR = wR(4); uR = wR(5); pR = wR(7); s11R = wR(9);
zR = rR*cR*cR + (4/3)*bR/(rR*rR);

% Initial guess
usg = ( (pL-s11L)*cL - (pR-s11R)*cR + zR*uR + zL*uL ) / (zR + zL);

% Iteration to find the u velocity state
us = fzero(@(us) ps_find(uL,uR,us,wL,wR),usg);
% us = usg;

% Pressure and s11 stress for left and right state for output
[psL,s11sL,rsL,s22sL] = ps_calc(uL,us,us,wL,-1);
[psR,s11sR,rsR,s22sR] = ps_calc(us,uR,us,wR,1);

%  storing the left state
wL(4) = rsL; wL(5) = us; wL(7) = psL; wL(9) = s11sL; wL(11) = s22sL;
% storing the right state
wR(4) = rsR; wR(5) = us; wR(7) = psR; wR(9) = s11sR; wR(11) = s22sR;

end

% function to find the zero such that the total stress is constant across
% contact
function [dsig] = ps_find(uL,uR,us,wL,wR)
    [psL,s11sL,~] = ps_calc(uL,us,us,wL,-1);
    [psR,s11sR,~] = ps_calc(us,uR,us,wR,1);
    sigR = psR-s11sR;
    sigL = psL-s11sL;
    dsig = sigL - sigR ;
end

% identifying if it is a shock or rarefaction
function [ps,s11s,rs,s22s] = ps_calc(uL,uR,us,w0,sign)
    if uL > uR
        [ps,s11s,rs,s22s] = ps_shock_calc(us,w0,sign);
    else
        [ps,s11s,rs,s22s] = ps_rare_calc(us,w0,sign);
    end
end

% shock conditions solver
function [ps,s11s,rs,s22s] = ps_shock_calc(us,w0,sign)
    g0 = w0(1);
    c0 = w0(2);
    b0 = w0(3);
    r0 = w0(4);
    u0 = w0(5);
    p0 = w0(7);
    pi0= w0(8);
    s110 = w0(9);
    s220 = w0(11);
    mu = (g0-1)/(g0+1);
    ap = 0.5*(g0+1)*r0*(us-u0)*(us-u0)/(p0+pi0);
    bp = (8*b0/3)/((g0+1)*r0*(p0+pi0));
    a = mu-2+(bp-ap);
    b = 1-2*mu-2*(mu*ap+bp);
    c = mu-mu*mu*ap+bp;
    phivec = roots([1 a b c]);
    phi = phivec((phivec>1));
    ps = phi*(p0+pi0) - pi0;
    m = sqrt( r0*r0*c0*c0*(((g0+1)/(2*g0))*phi+(g0-1)/(2*g0))+(4/3)*b0 );
    s11s = s110 + sign*(4/3)*(b0/m)*(u0-us);    
    rs = r0*((g0+1)*(ps+pi0)+(g0-1)*(p0+pi0))/((g0-1)*(ps+pi0)+(g0+1)*(p0+pi0));
    s22s = s220 - sign*(2/3)*(b0/m)*(u0-us);    
end

% rarefaction solver
function [ps,s11s,rs,s22s] = ps_rare_calc(us,w0,sign)
    g0 = w0(1);
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
    ifps = @(ps) (u0 - us) + sign*integral(fps,p0+pi0,ps+pi0);
    ps = fzero(ifps,p0);
    xi = ((ps+pi0)/(p0+pi0))^(-1/g0);
    s11s = s110 - (4*b0/(3*r0))*(1-xi);
    xi = ((ps+pi0)/(p0+pi0))^(1/g0);
    rs = r0*xi;
    s22s = s220 - (2*b0/(3*r0))*(1-xi);
end