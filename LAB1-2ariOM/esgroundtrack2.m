clear all
close all
clc
th_g0= 0;
om_E = deg2rad(15.04/3600);
mu_E = astroConstants(13);
%% case 1
n_orb =15;

a = 8350;
e = 0.19760;
i =deg2rad(60);
OM =deg2rad(270);
om =deg2rad(45);
th =deg2rad(230);
[r0,v0]=kep2car(a, e, i, OM, om, th, mu_E);
% original
[alpha,delta,lon,lat] = groundtrack( r0,v0,th_g0,n_orb)

%repeating
hold on;
k = 12; m = 1;
n = om_E*k/m;
a_rep=(mu_E/n^2).^(1/3);
[r00,v00]=kep2car(a_rep, e, i, OM, om, th, mu_E);
[alphamod,deltamod,lonmod,latmod] = groundtrackp( r00,v00,th_g0,n_orb,1)

%% case 2
n_orb =30;

a = 26600;
e = 0.74;
i =deg2rad(63.4);
OM =deg2rad(50);
om =deg2rad(280);
th =deg2rad(0);
[r0,v0]=kep2car(a, e, i, OM, om, th, mu_E);
% original
[alpha,delta,lon,lat] = groundtrack( r0,v0,th_g0,n_orb)

%repeating
hold on;
k = 2; m = 1;
n = om_E*k/m;
a_rep=(mu_E/n^2).^(1/3);
[r00,v00]=kep2car(a_rep, e, i, OM, om, th, mu_E);
[alphamod,deltamod,lonmod,latmod] = groundtrackp( r00,v00,th_g0,n_orb,1)

%% case 3.1
n_orb =30;

a = 7171.01;
e = 0;
i =deg2rad(0);
OM =deg2rad(0);
om =deg2rad(40);
th =deg2rad(0);
[r0,v0]=kep2car(a, e, i, OM, om, th, mu_E);
% original
[alpha,delta,lon,lat] = groundtrack( r0,v0,th_g0,n_orb)

%repeating
hold on;
k = 20; m = 2;
n = om_E*k/m;
a_rep=(mu_E/n^2).^(1/3);
[r00,v00]=kep2car(a_rep, e, i, OM, om, th, mu_E);
[alphamod,deltamod,lonmod,latmod] = groundtrackp( r00,v00,th_g0,n_orb,1)
%% case 3.2
i =deg2rad(30);

[r0,v0]=kep2car(a, e, i, OM, om, th, mu_E);
% original
[alpha,delta,lon,lat] = groundtrack( r0,v0,th_g0,n_orb)

%repeating
hold on;
k = 29; m = 2;
n = om_E*k/m;
a_rep=(mu_E/n^2).^(1/3);
[r00,v00]=kep2car(a_rep, e, i, OM, om, th, mu_E);
[alphamod,deltamod,lonmod,latmod] = groundtrackp( r00,v00,th_g0,n_orb,1)

%% case 3.3
i =deg2rad(98);

[r0,v0]=kep2car(a, e, i, OM, om, th, mu_E);
% original
[alpha,delta,lon,lat] = groundtrack( r0,v0,th_g0,n_orb)

%repeating
hold on;
k = 15; m = 1;
n = om_E*k/m;
a_rep=(mu_E/n^2).^(1/3);
[r00,v00]=kep2car(a_rep, e, i, OM, om, th, mu_E);
[alphamod,deltamod,lonmod,latmod] = groundtrackp( r00,v00,th_g0,n_orb,1)

%% case GEO
%% case 3.1
n_orb =500;

a_rep = 42166.167;
e = 0;
i =deg2rad(0);
OM =deg2rad(0);
om =deg2rad(0);
th =deg2rad(20);
[r00,v00]=kep2car(a_rep, e, i, OM, om, th, mu_E);
[alphamod,deltamod,lonmod,latmod] = groundtrackp( r00,v00,th_g0,n_orb,1)











