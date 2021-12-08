%% HW 4 Problem 1ii Code
%% Email questions to mrdz@umich.edu
clear all; close all; clc;
WL = [1, 0, 1];
WR = [0.125, 0, 0.1];
g = 1.4;
pguess = 0.5;
tol = 1e-15;
ps = root_find(pguess,tol,g,WL,WR);
LX0 = 0; LX1 = 1; dloc = 0.5; N = 7500; T = 0.2;
cfl = 0.95; NX = 200;
%plot exact solution
figure(1)
plot_rp(ps, g, WL, WR, LX0, LX1, dloc, N, T);
v1 = 'Exact';
L = legend(v1);
set(L,'Interpreter','Latex','FontSize',14);


