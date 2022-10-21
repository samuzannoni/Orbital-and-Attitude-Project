%% ES groundtrack
clear all
close all
clc
% NB: tutti i grafici,rispetto alle slide, sono posticipati start e end di
% pi per come è fatta la funzione grountrack
%% case 1
th_g0= 0;
n_orb =3.25;

r0 =[ -4578.219; -801.084; -7929.708];
v0 = [0.800; -6.037;1.385 ];

[alpha, delta, lon, lat] = groundtrack(r0,v0,th_g0,n_orb);

%% case 2
n_orb = 30;
r0 =[ 3108.128; -1040.299; -6090.022];
v0 = [5.743; 8.055; 1.555];

[alpha, delta, lon, lat] = groundtrack(r0,v0,th_g0,n_orb);

%% case 3.1
n_orb= 5;
r0 = [5493.312; 4609.436; 0.000]; 
v0 = [-4.792; 5.711; 0.000];

[alpha, delta, lon, lat] = groundtrack(r0,v0,th_g0,n_orb);

%% case 3.2
r0 = [5493.312; 3991.889; 2304.718 ];
v0 = [-4.792; 4.946; 2.856];

[alpha, delta, lon, lat] = groundtrack(r0,v0,th_g0,n_orb);

%% case 3.3
r0 = [5493.312; -641.510; 4564.578];
v0 = [-4.792; -0.795; 5.656];
  
[alpha, delta, lon, lat] = groundtrack(r0,v0,th_g0,n_orb);



