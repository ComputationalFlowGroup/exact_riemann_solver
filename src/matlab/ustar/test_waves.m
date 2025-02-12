%% 2015 JCP 5-WAVE RIEMANN SOLVER
% Email questions to mrdz@umich.edu
clear all; close all; clc;
% density times shear modulus, must be constant
beta = 1e13;
% left state
gL = 4.4;
rL = 1000;
uL = 0;
vL = 100;
pL = 1E3*1E5;
piL = 6e8; 
s11L = 0;
s12L = 0;
s22L = 0;
bL = beta;
cL = sqrt(gL*(pL+piL)/rL);
% right state
gR = 4.4;
rR = 1000;
uR = 0;
vR = -100;
pR = 1e5;
piR = 6e8;
s11R = 0;
s12R = 0;
s22R = 0;
bR = beta;
cR = sqrt(gR*(pR+piR)/rR);
% generating the left and right state vectors
wL = [gL,cL,bL,rL,uL,vL,pL,piL,s11L,s12L,s22L];
wR = [gR,cR,bR,rR,uR,vR,pR,piR,s11R,s12R,s22R];
% evaluating the star states
[wsL,wsR] = rp_sfind(wL,wR);
LX0 = 0; LX1 = 1; dloc = 0.5; N = 1000; T = 64e-6;
hold on;
rp_plot(wL, wsL, wR, wsR, LX0, LX1, dloc, N, T,...
    'f4riemann_5wavegav_sol.dat','cfd_d1nx200p4ao1am0rm3ga10.dat','.b')

%% GAV WAVE
clear all; close all; clc;
% density times shear modulus, must be constant
beta = 1e12;
% left state
gL = 4.4;
rL = 1000;
uL = 0;
vL = 100;
pL = 1E5*1E5;
piL = 6e8; 
s11L = 0;
s12L = 0;
s22L = 0;
bL = beta;
cL = sqrt(gL*(pL+piL)/rL);
% right state
gR = 4.4;
rR = 1000;
uR = 0;
vR = -100;
pR = 1e5;
piR = 6e8;
s11R = 0;
s12R = 0;
s22R = 0;
bR = beta;
cR = sqrt(gR*(pR+piR)/rR);
% generating the left and right state vectors
wL = [gL,cL,bL,rL,uL,vL,pL,piL,s11L,s12L,s22L];
wR = [gR,cR,bR,rR,uR,vR,pR,piR,s11R,s12R,s22R];
% evaluating the star states
[wsL,wsR] = rp_sfind(wL,wR);
LX0 = 0; LX1 = 1; dloc = 0.5; N = 1000; T = 64e-6;
hold on;
rp_plot(wL, wsL, wR, wsR, LX0, LX1, dloc, N, T,...
    'f4riemann_5wavegav_sol.dat','cfd_d1nx200p4ao1am0rm3ga9.dat','.b')

%% WEAK Impact
clear all; close all; clc;
% density times shear modulus, must be constant
beta = 1e13;
% left state
gL = 4.4;
rL = 1000;
uL = 10;
vL = 0*100;
pL = 1*1E5;
piL = 6e8; 
s11L = 0;
s12L = 0;
s22L = 0;
bL = beta;
cL = sqrt(gL*(pL+piL)/rL);
% right state
gR = 4.4;
rR = 1000;
uR = -10;
vR = -0*100;
pR = 1e5;
piR = 6e8;
s11R = 0;
s12R = 0;
s22R = 0;
bR = beta;
cR = sqrt(gR*(pR+piR)/rR);
% generating the left and right state vectors
wL = [gL,cL,bL,rL,uL,vL,pL,piL,s11L,s12L,s22L];
wR = [gR,cR,bR,rR,uR,vR,pR,piR,s11R,s12R,s22R];
% evaluating the star states
[wsL,wsR] = rp_sfind(wL,wR);
LX0 = 0; LX1 = 1; dloc = 0.5; N = 1000; T = 64e-6;
hold on;
rp_plot(wL, wsL, wR, wsR, LX0, LX1, dloc, N, T,...
    'f4riemann_5wavegav_sol.dat','cfd_d1nx200p4ao1am0rm3ga10.dat','.b')

%% Strong Impact
clear all; close all; clc;
% density times shear modulus, must be constant
beta = 1e13;
% left state
gL = 4.4;
rL = 1000;
uL = 1000;
vL = 0*100;
pL = 1*1E5;
piL = 6e8; 
s11L = 0;
s12L = 0;
s22L = 0;
bL = beta;
cL = sqrt(gL*(pL+piL)/rL);
% right state
gR = 4.4;
rR = 1000;
uR = -1000;
vR = -0*100;
pR = 1e5;
piR = 6e8;
s11R = 0;
s12R = 0;
s22R = 0;
bR = beta;
cR = sqrt(gR*(pR+piR)/rR);
% generating the left and right state vectors
wL = [gL,cL,bL,rL,uL,vL,pL,piL,s11L,s12L,s22L];
wR = [gR,cR,bR,rR,uR,vR,pR,piR,s11R,s12R,s22R];
% evaluating the star states
[wsL,wsR] = rp_sfind(wL,wR);
LX0 = 0; LX1 = 1; dloc = 0.5; N = 1000; T = 41e-6;
hold on;
rp_plot(wL, wsL, wR, wsR, LX0, LX1, dloc, N, T,...
    'f4riemann_5wavegav_sol.dat','cfd_d1nx200p4ao1am0rm3ga10.dat','.b')

%% Weak Exp
clear all; close all; clc;
% density times shear modulus, must be constant
beta = 1e13;
% left state
gL = 4.4;
rL = 1000;
uL = -50;
vL = 0*100;
pL = 1*1E5;
piL = 6e8; 
s11L = 0;
s12L = 0;
s22L = 0;
bL = beta;
cL = sqrt(gL*(pL+piL)/rL);
% right state
gR = 4.4;
rR = 1000;
uR = 50;
vR = -0*100;
pR = 1e5;
piR = 6e8;
s11R = 0;
s12R = 0;
s22R = 0;
bR = beta;
cR = sqrt(gR*(pR+piR)/rR);
% generating the left and right state vectors
wL = [gL,cL,bL,rL,uL,vL,pL,piL,s11L,s12L,s22L];
wR = [gR,cR,bR,rR,uR,vR,pR,piR,s11R,s12R,s22R];
% evaluating the star states
[wsL,wsR] = rp_sfind(wL,wR);
LX0 = 0; LX1 = 1; dloc = 0.5; N = 1000; T = 64e-6;
hold on;
rp_plot(wL, wsL, wR, wsR, LX0, LX1, dloc, N, T,...
    'f4riemann_5wavegav_sol.dat','cfd_d1nx200p4ao1am0rm3ga10.dat','.b')

%% Strong Exp
clear all; close all; clc;
% density times shear modulus, must be constant
beta = 1e13;
% left state
gL = 4.4;
rL = 1000;
uL = -800;
vL = 0*100;
pL = 1*1E5;
piL = 6e8; 
s11L = 0;
s12L = 0;
s22L = 0;
bL = beta;
cL = sqrt(gL*(pL+piL)/rL);
% right state
gR = 4.4;
rR = 1000;
uR = 800;
vR = -0*100;
pR = 1e5;
piR = 6e8;
s11R = 0;
s12R = 0;
s22R = 0;
bR = beta;
cR = sqrt(gR*(pR+piR)/rR);
% generating the left and right state vectors
wL = [gL,cL,bL,rL,uL,vL,pL,piL,s11L,s12L,s22L];
wR = [gR,cR,bR,rR,uR,vR,pR,piR,s11R,s12R,s22R];
% evaluating the star states
[wsL,wsR] = rp_sfind(wL,wR);
LX0 = 0; LX1 = 1; dloc = 0.5; N = 1000; T = 40e-6;
hold on;
rp_plot(wL, wsL, wR, wsR, LX0, LX1, dloc, N, T,...
    'f4riemann_5wavegav_sol.dat','cfd_d1nx200p4ao1am0rm3ga10.dat','.b')
