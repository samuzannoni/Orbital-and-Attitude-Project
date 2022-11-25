clc
clear all
close all

%% Task 1

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

T_sun = 365*24*3600; 
n_sun = 2*pi/T_sun;
epsilon = deg2rad(23.45); % Ecliptic's inclination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler Equation Solver

% Principal moments of inertia
Ix = 0.04; % Kgm^2 
Iy = 0.06; % Kgm^2
Iz = 0.08; % Kgm^2
I = [Ix 0 0; 0 Iy 0; 0 0 Iz]; % Inertia matrix
I_inv = inv(I);

v1 = [0 1 -1]';
v2 = [-1 0 1]';
v3 = [1 -1 0]';

A0 =  eye(3); % initial DCM
wx_0 = 1e-06; % [Rad/s]
wy_0 = 1e-06; % [Rad/s]
wz_0 = n; % Mean motion of spacecraft, if satellite always gives the same face to Earth, w_z = n [Rad/s]
w_0 = [wx_0 wy_0 wz_0]'; % initial spin velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Task 2 and 3

% Inertia matrix of deployed elements (solar arrays)
J_depl = 1e-2*[100.9 0 0; 0 25.1 0; 0 0 91.6]; % [kg*m^2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacecraft main body data

% size of spacecraft main body
base = 0.263;
height = 0.366;
% normal to surfaces
n1 = [1;0;0];
n2 = [0;1;0];
n3 = [0;0;1];
n4 = -n1;
n5 = -n2;
n6 = -n3;
% surfaces
A1 = base*height;
r1 = [base/2; 0; 0];
A2 = A1;
A3 = base*base;
A4 = A1;
A5 = A2;
A6 = A3;
% distance from CG to centre of pressure on each surface
r_mb = [base/2 0 0; 0 base/2 0; 0 0 height/2; -base/2 0 0; 0 -base/2 0; 0 0 -height/2];
% radiation coefficients
rho_s_mb = 0.5; % specular reflection coeff of main body
rho_d_mb = 0.1; % diffuse reflection coeff of main body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solar panels data (two panels)

% Panel 1 data:
% size of solar panel #1
height_p1 = 0.1;
length_p1 = 0.7;
% normal to surfaces
n7=[1;0;0];
n8 = -n7;
% surfaces
A7 = length_p1*height_p1;
A8 = A7;

% Panel 2 data:
% size of solar panel #2
length_p2 = length_p1; % size of solar panel #2
height_p2 = height_p1;
% normal to surfaces
n9 = [1;0;0];
n10 = -n9;
% surfaces
A9 = length_p1*height_p1;
A10 = A9;

% distance from CG to centre of pressure on each surface
r_sp = [0 length_p1/2+base/2 0; 0 length_p1/2+base/2 0; 0 -length_p2/2-base/2 0; 0 -length_p2/2-base/2 0];
% radiation coefficients
rho_s_sp = 0.1; % specular reflection of solar arrays
rho_d_sp = 0.1; % diffuse reflection of solar arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fe = 1358; % W/m^2
c = 3*1e8; % [m/s], speed of light
P = Fe/c; % Solar radiation pressure

Area_mb = diag([A1 A2 A3 A4 A5 A6]);
Area_sp = diag([A7 A8 A9 A10]);
N_b_mb  = [n1'; n2'; n3'; n4'; n5'; n6']; % matrix of vectors n for main body
N_b_sp  = [n7'; n8'; n9'; n10']; % matrix of vectors n for solar panels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Task 4

omega_LN = [0 0 n]'; 
unit = [1 0 0]'; % vector_C = A_B/L*unit (see slides)
i_matrix = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]; % matrix that multiply A_B/L to get the LVLH with an inclined orbit (see slides)

% WARNING: must add inertia of the deployed elements to the inertia
% matrixes!
%% Task 5

% Coefficients according to IGRF2000
g1_0 = -29615*1e-9;                     %[T]
g1_1 = -1728*1e-9;                      %[T]
h1_1 = 5186*1e-9;                       %[T]
H0 = ( g1_0^2 + g1_1^2 + h1_1^2 )^(1/2);

m = [0.01; 0.05; 0.01]; % parasitic magnetic torque vector: considered a small spacecraft, maximum m (worst case) is [0.1 0.1 0.1]'

%% simulation

simulation = sim("Lab07.slx");
T_MAG = simulation.simout;
T_GG = simulation.simout1;
T_SRP = simulation.simout2;
c = simulation.c_vector;
A = simulation.simout3;