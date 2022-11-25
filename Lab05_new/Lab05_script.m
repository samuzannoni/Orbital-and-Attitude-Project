clc
clear all
close all

%% Pointing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
mu_E = astroConstants(13); % Gravitational constant Earth
R_E = astroConstants(23); % Earth's radius
Rs = astroConstants(2); % Distance Earth-Sun
alpha = deg2rad(11.5); % Inclination of Magnetic Axis
w_E = deg2rad(15.04)/3600; % Earth's rotational speed, [rad/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Model

a = 16733.65;
e = 0.19760;
i = deg2rad(60);
om = deg2rad(45);
theta0 = deg2rad(230);
n = sqrt(mu_E/(a^3)); % Mean motion satellite
T = 2*pi/n; % Period of orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler Equation Solver

% Principal moments of inertia
Ix = 0.07; % Kgm^2 
Iy = 0.0504; % Kgm^2
Iz = 0.0109; % Kgm^2
I = [Ix 0 0; 0 Iy 0; 0 0 Iz]; % Inertia matrix
I_inv = inv(I);

v1 = [0 1 -1]';
v2 = [-1 0 1]';
v3 = [1 -1 0]';

A0 =  eye(3); % initial DCM
wx_0 = 0.5; % [Rad/s]
wy_0 = 0.1; % [Rad/s]
wz_0 = 0.1; % [Rad/s]
w_0 = [wx_0; wy_0; wz_0]; % initial spin velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target pointing performance

h_0_N = [1; 0; 0]; % initial attitude parameters in inertial frame
h_N_target = h_0_N; % satellite must point z-axis in inertial frame at each timestep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of attitude parameters at initial timestep

% oss: attitude parameters depend upon angular velocity so i must start
% from h_0_target then compute A_B/N target at t=t0

A_bn_target_t0_transposed = ((I'.*w_0').*h_N_target)./((norm(h_N_target))^2);

%% simulation

simulation = sim("Lab05.slx");
h_N = simulation.h_inertial;
A_BN = simulation.matrix;
err = simulation.simout;

