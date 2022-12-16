clc
clear all
close all
format long g;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

mu_E = astroConstants(13); % Gravitational constant Earth
R_E = astroConstants(23); % Earth's radius
Rs = astroConstants(2); % Distance Earth-Sun
alpha = deg2rad(11.5); % Inclination of Magnetic Axis
w_E = deg2rad(15.04)/3600; % Earth's rotational speed, [rad/s]
epsilon = deg2rad(23.45); % Ecliptic inclination
T_sun = 365*24*3600; 
n_sun = 2*pi/T_sun;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbital Model: circular LEO (STARLINK-1053)

% TLE:
% 1 44758C 19074AX  22350.05604167 -.00020847  00000-0 -13966-2 0  3501
% 2 44758  53.0528 122.2873 0001456  77.1723  38.0377 15.06422412    17

a = R_E + 540; % orbital semi-major axis [km]
e = 0; % orbital eccentricity
i = deg2rad(86.4); % orbital inclination
om = deg2rad(77.1723); % argument of pericenter
theta0 = deg2rad(0); % true anomaly
n = sqrt(mu_E/(a^3)); % angular velocity
T = 2*pi/n; % orbital period

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Euler Equation Solver

% Principal moments of inertia
Ix = 0.36; % Kgm^2 
Iy = 0.36; % Kgm^2
Iz = 0.01; % Kgm^2
I = [Ix 0 0; 0 Iy 0; 0 0 Iz]; % Inertia matrix
I_inv = inv(I);

v1 = [0 1 -1]';
v2 = [-1 0 1]';
v3 = [1 -1 0]';

% Initial conditions:

A0 =  eye(3); % initial DCM

wx_0 = 1e-06; % [Rad/s]
wy_0 = 1e-06; % [Rad/s]
wz_0 = n; % Mean motion of spacecraft, if satellite always gives the same face to Earth, w_z = n [Rad/s]
w_0 = [wx_0 wy_0 wz_0]';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spacecraft main body data: Cubesat 12U 120x10x10 cm, 1 kg

% Obs.: mass needs to be modified after the choice of sensors and actuators!
% Assumptions: Consider solar panels attached on spacecraft body!
%               x towards Earth, y towards clockwise motion (opposite to motion of spacecraft), z towards North

% size of spacecraft main body
base = 0.1; % [m]
height = 12*base; % [m]
% normal to surfaces
n1 = [1;0;0];
n2 = [0;1;0];
n3 = [0;0;1];
n4 = -n1;
n5 = -n2;
n6 = -n3;
% surfaces
A1 = base*height;
A2 = A1;
A3 = base*base;
A4 = A1;
A5 = A2;
A6 = A3;
% distance from CG to centre of pressure on each surface
r_mb = [base/2 0 0; 0 base/2 0; 0 0 height/2; -base/2 0 0; 0 -base/2 0; 0 0 -height/2];
% radiation coefficients
rho_s_mb = 0.1; % specular reflection coeff of main body
rho_d_mb = 0.1; % diffuse reflection coeff of main body

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deployed elements data

% % Inertia matrix of deployed elements (solar arrays)
% J_depl = 1e-2*[100.9 0 0; 0 25.1 0; 0 0 91.6]; % [kg*m^2]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Solar panels data (two panels)
% 
% % Panel 1 data:
% % size of solar panel #1
% height_p1 = 0.1;
% length_p1 = 0.7;
% % normal to surfaces
% n7=[1;0;0];
% n8 = -n7;
% % surfaces
% A7 = length_p1*height_p1;
% A8 = A7;
% 
% % Panel 2 data:
% % size of solar panel #2
% length_p2 = length_p1; % size of solar panel #2
% height_p2 = height_p1;
% % normal to surfaces
% n9 = [1;0;0];
% n10 = -n9;
% % surfaces
% A9 = length_p1*height_p1;
% A10 = A9;
% 
% % distance from CG to centre of pressure on each surface
% r_sp = [0 length_p1/2+base/2 0; 0 length_p1/2+base/2 0; 0 -length_p2/2-base/2 0; 0 -length_p2/2-base/2 0];
% % radiation coefficients
% rho_s_sp = 0.1; % specular reflection of solar arrays
% rho_d_sp = 0.1; % diffuse reflection of solar arrays

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solar Radiation Pressure Torque --> it's negligible! order: o(10^-25)

Fe = 1358; % W/m^2
c = 3*1e8; % [m/s], speed of light
P = Fe/c; % Solar radiation pressure

Area_mb = diag([A1 A2 A3 A4 A5 A6]);
% Area_sp = diag([A7 A8 A9 A10]);
N_b_mb  = [n1'; n2'; n3'; n4'; n5'; n6']; % matrix of vectors n for main body
% N_b_sp  = [n7'; n8'; n9'; n10']; % matrix of vectors n for solar panels

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gravity Torque: o(10^-7)

omega_LN = [0 0 n]'; 
unit = [1 0 0]'; % vector_C = A_B/L*unit (see slides)
i_matrix = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]; % matrix that multiply A_B/L to get the LVLH with an inclined orbit (see slides)

% When spacecrafgt spins around z-axis, it must remain stable even with torque applied
% because z-axis is major or minor inertial axis (single spin stability). Use format long. 

%% Magnetic Torque: o(10^-7)

% Coefficients according to IGRF2000
g1_0 = -29615*1e-9;                     %[T]
g1_1 = -1728*1e-9;                      %[T]
h1_1 = 5186*1e-9;                       %[T]
H0 = ( g1_0^2 + g1_1^2 + h1_1^2 )^(1/2);

m = [0.01; 0.05; 0.01]; % parasitic magnetic torque vector: considered a small spacecraft, maximum m (worst case) is [0.1 0.1 0.1]'


%% Star sensor (mounted on face A4, n4 = [-1; 0; 0]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensor Data:
% Sensor chosen is "Sagitta Star Tracker" (it's a star mapper)
% FOV = 23.4 deg
% Limit tracked magnitude: [0:1]
% Accuracy:
    cross_boresight = 2; % ["]
    around_boresight = 10;  % ["]

    cross_boresight = deg2rad(dms2degrees([0,0,cross_boresight])); % [rad]
    around_boresight = deg2rad(dms2degrees([0,0,around_boresight])); % [rad]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumptions:
% GEO orbit (i=0) above Kourou: 0 deg, 0' 0.0''; W 52deg 39' 01"
% Day in which the attitude is analysed: 2 Jan 2023

% I assume located stars are very near to the local meridian of the
% satellite, so we can assume LST = RA (local sidereal time equal to right
% ascension of the star)...this assumption could be negligible!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Catalogue of stars (ref J2000 epoch, sky seen on 2 Jan 2023 from equator above Kourou (0 deg, 0' 0.0''; W 52deg 39' 01'')):
%
%       Comment: take 6 stars spread along different regions of the sky so as to cover the whole orbit
%
% Altair: 19h 50m 46.51s; 8deg 52' 14.3''; Mag: 0.75; 16.73 ly 
%       WARNING: Altair is quite close to the Sun!
% Spica: 13h 25m 11.57s; -11deg 09' 40.75''; Mag: 0.95; 250 ly
% Procyon:  7h 39m 18,1183s; 5deg 13′ 29,975″; Mag: 0.40; 11.4 ly
% Rigel: 05h 14m 32,30s; -08deg 12′ 06,89″; Mag: 0.15; 860 ly
% Betelgeuse: 05h 55m 10,3053s; 07deg 24′ 25,426″; Mag: 0.45; 500 ly
% Aldebaran: 4h 35m 55,239s; 16deg 30′ 33,488″; Mag: 0.85; 66.64 ly
%       Comment: No bright star between R.A. 20h and 4h along the equator
% Discarded stars (too low magnitude):
% Vega: 18h 36m 56,3364s; 38deg 47′ 1,291″; Mag 0.03; 25.3 ly
% Capella: 5h 16m 41,359s; +45deg 59′ 52,768″; Mag: 0.08; 43 ly
% Rigil Kentaurus: 14h 39m 36,495s; −60deg 50′ 2,308″; Mag: -0.27; 4.364 ly
% Canopus: 06h 23m 57.10988s; 52deg 41′ 44.3810″; Mag: -0.74; 320 ly
% Arcturus: 14h 15m 39.7s; 19deg 10′ 56″; Mag: -0.05; 36.2 ly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the catalogue:
% R.A in first column, DEC in second column, distance [km] in third column

n_stars = 6;
catalogue=zeros(n_stars,3);
catalogue(:,1) = deg2rad([time2deg(19,50,46.51); time2deg(13,25,11.57); time2deg(7,39,18.1183); time2deg(5,14,32.30); time2deg(5,55,10.3053); time2deg(4,35,55.239)]); % [rad]
catalogue(:,2) = deg2rad(dms2degrees([8,52,14.3; -11,9,40.75; 5,13,29.975; -8,12,6.89; 7,24,25.426; 16,30,33.488])); % [rad]
catalogue(:,3) = [ly2km(16.73); ly2km(250); ly2km(11.4); ly2km(860); ly2km(500); ly2km(66.64)]; % [km]

m_N = zeros(size(catalogue,1),size(catalogue,2));
for i=1:size(catalogue, 1)

    [m_N(i,1),m_N(i,2),m_N(i,3)] = sph2cart(catalogue(i,1),catalogue(i,2),catalogue(i,3)); % WARNING: In function sph2cart angles must be in radians

end

m_N = m_N'; % I transpose m_N to get the location of each star on every column ( m_N: [n_starsx3] -> m_N':[3xn_stars] )
magnitude_columns_m_N = vecnorm(m_N);

m_N_unit = zeros(size(m_N,1),size(m_N,2));
for j=1:length(magnitude_columns_m_N)
    
    m_N_unit(:,j) = m_N(:,j)./magnitude_columns_m_N(j);
end

% I must multiply A_B/N per a unit vector: that's why I transfromed m_N into a unit vector (m_N_unit) !

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error matrix is built in Simulink:

%       Random number block parameters:
%           Mean: 0
%           Variance: 0.16
%           Seed: 0
%           Sample time: 0.1

% The following code for the genration of random numbers was found on: https://it.mathworks.com/help/matlab/math/random-numbers-with-specific-mean-and-variance.html
% rng('shuffle','twister');
% a = 0.4; % variance is a^2: 0.4^2 = 0.16
% b = 0; % mean value
% random_error = a.*randn + b; % generate a random error factor with normal distribution with mean 0 and variance 0.16 so as to belong to [-1;+1]

% A_epsilon = [1, around_boresight/2 + around_boresight/2*time_effect, -(cross_boresight/2 + cross_boresight/2*time_effect); ...
%              -(around_boresight/2 + around_boresight/2*time_effect), 1, cross_boresight/2 + cross_boresight/2*time_effect; ...
%              cross_boresight/2 + cross_boresight/2*time_effect, -(cross_boresight/2 + cross_boresight/2*time_effect), 1]; % [rad]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Attitude determination

peso = 1/6;

h_N_target = [0; 0; 1]; % Target angular momentum (it points North)

%% Actuators

%% Control



%% simulation

simulation = sim("Spacecraft_Model.slx");
% m_B = simulation.simout;
% torque_srp = simulation.srp;
% torque_gg = simulation.gg;
% torque_mag = simulation.mag;