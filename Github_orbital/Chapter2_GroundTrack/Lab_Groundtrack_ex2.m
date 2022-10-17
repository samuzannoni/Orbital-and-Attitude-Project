clc
clear all
close all

% First orbit
r0 = [-4578.219, -801.084, -7929.708]'; % Km
v0 = [0.800, -6.037, 1.385]'; % Km/s

a = 8350; % Km
e = 0.1976; % eccentricity
i = 60; % deg 
Om = 270; % deg 
om = 45; % deg 
f0 = 230; % deg: initial true anomaly

mu_E = astroConstants(13);
T = 2*pi *sqrt(a^3/mu_E);
w_E=deg2rad(15.04/3600);
T_E=2*pi/w_E;
mk_ratio=T/T_E