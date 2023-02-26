% main script to the piece-wise constant Riemann problem for the 
% Euler equations exactly
clear all; close all; clc; % clearing all variables for work

% Consider that we are given a WL and WR primitive vectors
% all values are non-dimensional
% gamma
g = 1.4; 
tolerance = 1E-5; % this is the tolerance at which your root solver will stop
LX0 = 0; LX1 = 1; N = 1000; T = 0.25;
% case 1
% left state
WL = [1, 0, 1]; % first entry is the density, second the velocity, third the pressure
% right state
WR = [0.125, 0, 0.1]; % same as the left state WL
% initial pressure guess
pguess = 0.5*(WL(3)+WR(3)); % non-dimensional pressure
[pstar] = euler_riemann_root_finder(pguess,WL,WR,g,tolerance);
plot_riemann(pstar, g, WL, WR, LX0, LX1, N, T,1);
saveas(gca,'case1.png')

% case 2
% left state
WL = [1, -2, 0.4]; % first entry is the density, second the velocity, third the pressure
% right state
WR = [1, 2, 0.4]; % same as the left state WL
T = 0.15;
pguess = 0.5*(WL(3)+WR(3));
[pstar] = euler_riemann_root_finder(pguess,WL,WR,g,tolerance);
plot_riemann(pstar, g, WL, WR, LX0, LX1, N, T,2);
saveas(gca,'case2.png')

% case 3
% left state
WL = [1, 0, 1000]; % first entry is the density, second the velocity, third the pressure
% right state
WR = [1, 0, 0.01]; % same as the left state WL
T = 0.012;
pguess = 0.5*(WL(3)+WR(3));
[pstar] = euler_riemann_root_finder(pguess,WL,WR,g,tolerance);
plot_riemann(pstar, g, WL, WR, LX0, LX1, N, T,3);
saveas(gca,'case3.png')

% case 4
% left state
WL = [1, 0, 0.01]; % first entry is the density, second the velocity, third the pressure
% right state
WR = [1, 0, 100]; % same as the left state WL
T = 0.035;
pguess = 0.5*(WL(3)+WR(3));
[pstar] = euler_riemann_root_finder(pguess,WL,WR,g,tolerance);
plot_riemann(pstar, g, WL, WR, LX0, LX1, N, T,4);
saveas(gca,'case4.png')