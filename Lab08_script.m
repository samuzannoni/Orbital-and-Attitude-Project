clc
clear all
close all

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
%% Orbital Model: GEO (circular equatorial orbit)

a = 42164; % orgital semi-major axis [km]
e = 0; % orbital eccentricity
i = deg2rad(0); % orbital inclination
om = deg2rad(0); % argument of pericenter
theta0 = deg2rad(0); % true anomaly
n = sqrt(mu_E/(a^3)); % angular velocity
T = 2*pi/n; % orbital period

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Euler Equation Solver

% Principal moments of inertia
Ix = 0.04; % Kgm^2 
Iy = 0.06; % Kgm^2
Iz = 0.08; % Kgm^2
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
%% Spacecraft main body data:

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deployed elements data

% % Inertia matrix of deployed elements (solar arrays)
% J_depl = 1e-2*[100.9 0 0; 0 25.1 0; 0 0 91.6]; % [kg*m^2]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solar panels data (two panels)

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solar Radiation Pressure Torque

Fe = 1358; % W/m^2
c = 3*1e8; % [m/s], speed of light
P = Fe/c; % Solar radiation pressure

Area_mb = diag([A1 A2 A3 A4 A5 A6]);
Area_sp = diag([A7 A8 A9 A10]);
N_b_mb  = [n1'; n2'; n3'; n4'; n5'; n6']; % matrix of vectors n for main body
N_b_sp  = [n7'; n8'; n9'; n10']; % matrix of vectors n for solar panels

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gravity Torque

omega_LN = [0 0 n]'; 
unit = [1 0 0]'; % vector_C = A_B/L*unit (see slides)
i_matrix = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]; % matrix that multiply A_B/L to get the LVLH with an inclined orbit (see slides)

% When spacecrafgt spins around z-axis, it must remain stable even with torque applied
% because z-axis is major or minor inertial axis (single spin stability). Use format long. 

%% Magnetic Torque

% Coefficients according to IGRF2000
g1_0 = -29615*1e-9;                     %[T]
g1_1 = -1728*1e-9;                      %[T]
h1_1 = 5186*1e-9;                       %[T]
H0 = ( g1_0^2 + g1_1^2 + h1_1^2 )^(1/2);

m = [0.01; 0.05; 0.01]; % parasitic magnetic torque vector: considered a small spacecraft, maximum m (worst case) is [0.1 0.1 0.1]'


%% Star sensor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensor Data:
% Sensor chosen is "Sagitta Star Tracker" (it's a star mapper)
% FOV = 23.4 deg
% Limit tracked magnitude: [0:1]
% Accuracy:
    cross_boresight = 2; % ["]
    around_boresight = 10;  % ["]

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
% Spica: 13h 25m 11.57s; -11deg 09' 40.75''; Mag: 0.95; 250 ly
% Procyon:  7h 39m 18,1183s; 5deg 13′ 29,975″; Mag: 0.40; 11.4 ly
% Rigel: 05h 14m 32,30s; -08deg 12′ 06,89″; Mag: 0.15; 860 ly
% Betelgeuse: 05h 55m 10,3053s; 07deg 24′ 25,426″; Mag: 0.45; 500 ly
% Aldebaran: 4h 35m 55,239s; 16deg 30′ 33,488″; Mag: 0.85; 66.64 ly

% Discarded stars:
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
catalogue(:,1) = [time2deg(19,50,46.51); time2deg(13,25,11.57); time2deg(7,39,18.1183); time2deg(5,14,32.30); time2deg(5,55,10.3053); time2deg(4,35,55.239)];
catalogue(:,2) = dms2degrees([8,52,14.3; -11,9,40.75; 5,13,29.975; -8,12,6.89; 7,24,25.426; 16,30,33.488]);
catalogue(:,3) = [ly2km(16.73); ly2km(250); ly2km(11.4); ly2km(860); ly2km(500); ly2km(66.64)];

m_N = zeros(size(catalogue,1),size(catalogue,2));
for i=1:size(catalogue, 1)

    [m_N(i,1),m_N(i,2),m_N(i,3)] = sph2cart(catalogue(i,1),catalogue(i,2),catalogue(i,3));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the error matrix:

cross_boresight = deg2rad(dms2degrees([0,0,cross_boresight])); 
around_boresight = deg2rad(dms2degrees([0,0,around_boresight]));

time_effect = deg2rad(rand);

A_epsilon = [1, around_boresight/2 + around_boresight/2*time_effect, -(cross_boresight/2 + cross_boresight/2*time_effect); ...
             -(around_boresight/2 + around_boresight/2*time_effect), 1, cross_boresight/2 + cross_boresight/2*time_effect; ...
             cross_boresight/2 + cross_boresight/2*time_effect, -(cross_boresight/2 + cross_boresight/2*time_effect), 1]; % [rad]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% simulation
peso= 1/6;
V = m_N';
V_pseudoinv = V'*(inv(V*V'));

simulation = sim("Lab08.slx");
m_B = simulation.simout;

ABN = simulation.simout7;
errore= simulation.simout4;
% out1 = simulation.simout5;
% out2 = simulation.simout6;