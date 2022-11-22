clc
clear all
close all

%% Task 1

% Principal moments of inertia: case 1 
Ix = 0.04; % Kgm^2 
Iy = 0.06; % Kgm^2
Iz = 0.08; % Kgm^2
I = [Ix 0 0; 0 Iy 0; 0 0 Iz]; % Inertia matrix
I_inv = inv(I); % Matrice inversa di I

mu_E = astroConstants(13); % Gravitational constant Earth
R_E = astroConstants(23); %Earth's radius

% Orbit
a = 16733.65; % Km
e = 0.19760;
i = deg2rad(60);
om = deg2rad(45);
theta0 = deg2rad(230);
n = sqrt(mu_E/(a^3)); % mean motion
T = 2*pi/n; % Period

v1 = [0 1 -1]';
v2 = [-1 0 1]';
v3 = [1 -1 0]';

A0 =  eye(3);

wx_0 = 1e-06; % Rad/s
wy_0 = 1e-06; % Rad/s
wz_0 = n; % Rad/s

w_0 = [wx_0 wy_0 wz_0]'; % vettore condizioni iniziali

w_L = [0 0 n]';

T_sun = 365*24*3600; 
n_sun = 2*pi/T_sun;
epsilon = deg2rad(23.45); %ecliptic
Rs = astroConstants(2);

% simulation = sim("Lab7_task1");
% 
% theta = rad2deg(simulation.theta);
% t = simulation.tout;
% A_body = simulation.A_body;

%% Task 2

% inertia matrix of deployed elements (solar arrays)
J_depl = 1e-2*[100.9 0 0; 0 25.1 0; 0 0 91.6]; % [kg*m^2]

% Spacecraft main body
base = 0.263;
height = 0.366;

n1 = [1;0;0];
n2 = [0;1;0];
n3 = [0;0;1];
n4 = -n1;
n5 = -n2;
n6 = -n3;

% surface [m^2]
A1 = base*height;
A2 = A1;
A3 = base*base;
A4 = A1;
A5 = A2;
A6 = A3;

rho_s_mb = 0.5; % specular reflection of main body
rho_d_mb = 0.1; % diffuse reflection of main body

% Solar panels (two panels)
height_p1 = 0.1; % size of solar panel #1
length_p1 = 0.7;
n7=[1;0;0];
n8 = -n7;
A7 = length_p1*height_p1;
A8 = A7;

length_p2 = length_p1; % size of solar panel #2
height_p2 = height_p1;
n9 = [1;0;0];
A9 = length_p1*height_p1;
n10 = -n9;
A10 = A9;

rho_s_sp = 0.1; % specular reflection of solar arrays
rho_d_sp = 0.1; % diffuse reflection of solar arrays
%% Task 3

Fe = 1358; % W/m^2
c = 3*1e8; % [m/s], speed of light
P = Fe/c; % Solar radiation pressure

%Area_mb = diag([A1 A2 A3 A4 A5 A6]);
%Area_sp = diag([A7 A8 A9 A10]);
%N_b_mb  = [n1'; n2'; n3'; n4'; n5'; n6']; % matrix of vectors n for main body
%N_b_sp  = [n7'; n8'; n9'; n10']; % matrix of vectors n for solar panels

r1 = [base/2, 0, 0 ];
r2 = [ 0, base/2, 0];
r3 = [ 0, 0, height/2];
r4 = -r1;
r5 = -r2;
r6 = -r3;


% distance from CG to centre of pressure on each surface
r7 = [0 length_p1/2+base/2 0];
r8=[0 length_p1/2+base/2 0];
r9=[0 -length_p2/2-base/2 0];
r10=[0 -length_p2/2-base/2 0];
% radiation coefficients
simulation = sim("Lab7scalar");
% simulation = sim("Lab7");
% 
% S_B = simulation.S_B;
% cos_theta = simulation.cos_theta;
% ctrl = simulation.bucchina;
% F = simulation.Force;
M = simulation.Torque;



% try torque on solar arrays (1e-3<T<1e-6), torque on satellite faces could be very small
% rename the simulink model