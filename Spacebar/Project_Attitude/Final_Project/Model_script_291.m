clc
clear all
close all
format short g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

mu_E = astroConstants(13); % Gravitational constant Earth
R_E = astroConstants(23); % Earth's radius
Rs = astroConstants(2); % Distance Earth-Sun
alpha = deg2rad(11.5); % Inclination of Magnetic Axis
w_E = deg2rad(15.04)/3600; % Earth's rotational speed, [rad/s]
epsilon = deg2rad(23.45); % Ecliptic inclination
T_sun = 365*24*3600; % Revolution period of Earth (1 year) [s]
n_sun = 2*pi/T_sun; % Mean motion of Earth with respect to the Sun

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbital Model: circular LEO (STARLINK-1053)

% TLE:
%   1 44758C 19074AX  23003.86506944 -.00005627  00000-0 -37706-3 0    33
%   2 44758  53.0545  37.8515 0001824  45.9004 247.5767 15.06399717    15
% SOURCE: https://orbit.ing-now.com/satellite/44758/2019-074ax/starlink-1053/
%         https://www.n2yo.com/satellite/?s=44758
%         https://www.nanosats.eu/

a = R_E + 550; % Orbital semi-major axis [km]
e = 0.0001824; % Orbital eccentricity
i = deg2rad(53.0545); % Orbital inclination
om = deg2rad(45.9004); % Argument of pericenter
theta0 = deg2rad(0); % True anomaly
n = sqrt(mu_E/(a^3)); % Angular velocity
T = 2*pi/n; % Orbital period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spacecraft main body data: Cubesat 12U 120x10x10 cm, 1 kg

% Obs.: mass needs to be modified after the choice of sensors and actuators!
% Reference system: x towards Earth, y towards clockwise motion (opposite to motion of spacecraft), z towards North

% Size of spacecraft main body:
base = 0.1; % [m]
height = 12*base; % 12U [m]
mass_mb = 12; % [kg]

% Principal moments of inertia of spacecraft main body:
Ix_mb = mass_mb*(height/2)^2; % [Kg*m^2] 
Iy_mb = mass_mb*(height/2)^2; % [Kg*m^2]
Iz_mb = mass_mb*(base/2)^2; % [Kg*m^2]
I_mb = [Ix_mb 0 0; 0 Iy_mb 0; 0 0 Iz_mb]; % [Kg*m^2]

% Normal to surfaces:
n1 = [1;0;0];
n2 = [0;1;0];
n3 = [0;0;1];
n4 = -n1;
n5 = -n2;
n6 = -n3;

% Surfaces:
A1 = base*height;
A2 = A1;
A3 = base*base;
A4 = A1;
A5 = A2;
A6 = A3;

% Distance from CG to centre of pressure on each surface:
r_mb = [base/2 0 0; 0 base/2 0; 0 0 height/2; -base/2 0 0; 0 -base/2 0; 0 0 -height/2];

% Radiation coefficients:
rho_s_mb = 0.3; % specular reflection coeff of main body
rho_d_mb = 0.6; % diffuse reflection coeff of main body

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deployed elements data: two solar panels

% Install two panels along the y-direction on surface A2 and A5, with the dimensions given by slide 3 of "Labb07.pdf"

% Panel #1 data:

% size of solar panel #1
height_p1 = 1; % [m]
length_p1 = 5; % [m]

% SOURCE: https://en.wikipedia.org/wiki/Starlink

% normal to surfaces of solar panel #1
n7 = [1;0;0];
n8 = -n7;

% surfaces of solar panel #1
A7 = length_p1*height_p1;
A8 = A7;


% Panel #2 data:

% size of solar panel #2
length_p2 = length_p1; 
height_p2 = height_p1;

% normal to surfaces of solar panel #2
n9 = [1;0;0];
n10 = -n9;

% surfaces of solar panel #2
A9 = length_p1*height_p1;
A10 = A9;

% distance from CG to centre of pressure on each surface
r_sp = [0 length_p1/2+base/2 0; 0 length_p1/2+base/2 0; 0 -length_p2/2-base/2 0; 0 -length_p2/2-base/2 0];

% radiation coefficients
rho_s_sp = 0; % specular reflection of solar arrays
rho_d_sp = 0; % diffuse reflection of solar arrays

% Principal moments of inertia of spacecraft solar panels:
number_sp = 2;
mass_unit_surface_sp = 1.76; % [kg/(m^2)] % SOURCE: https://www.spectrolab.com/DataSheets/Panel/panels.pdf
surface_sp = height_p1*length_p1;
mass_sp = number_sp*(mass_unit_surface_sp*(surface_sp));

Ix_sp = mass_sp*(length_p1)^2; % [Kg*m^2] 
Iy_sp = mass_sp*(height_p1/2)^2; % [Kg*m^2]
Iz_sp = mass_sp*(length_p1)^2; % [Kg*m^2]
I_depl = [Ix_sp 0 0; 0 Iy_sp 0; 0 0 Iz_sp]; % [Kg*m^2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Euler Equation Solver

% Total inertia matrix (spacecraft main body + deployable elements):
I = I_mb + I_depl; % [kg*m^2]
I_inv = inv(I);

% Coefficients needed to solve Euler's equation:
v1 = [0 1 -1]';
v2 = [-1 0 1]';
v3 = [1 -1 0]';

% Initial conditions:
A0 = eye(3); % initial DCM

wx_0 = 0.9; % [Rad/s]
wy_0 = 0.7; % [Rad/s]
wz_0 = 0.3; % [Rad/s]

w_0 = [wx_0 wy_0 wz_0]'; % [Rad/s], initial velocity

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solar Radiation Pressure Torque --> it's negligible! order: o(10^-25)

% NOTE: It should be even lower because sspacecraft has a large eclipse time.

Fe = 1358; % W/m^2
c = 3*1e8; % [m/s], speed of light
P = Fe/c; % Solar radiation pressure

Area_mb = diag([A1 A2 A3 A4 A5 A6]);
Area_sp = diag([A7 A8 A9 A10]);
N_b_mb  = [n1'; n2'; n3'; n4'; n5'; n6']; % matrix of normal vectors for main body
N_b_sp  = [n7'; n8'; n9'; n10']; % matrix of normal vectors for solar panels

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gravity Torque: o(10^-4)

omega_LN = [0 0 n]'; 
unit = [1 0 0]'; % vector_C = A_B/L*unit (see slides)
i_matrix = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]; % matrix that multiply A_B/L to get the LVLH with an inclined orbit (see slides)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Magnetic Torque: o(10^-6)

% Coefficients according to IGRF2000
g1_0 = -29615*1e-9;                     %[T]
g1_1 = -1728*1e-9;                      %[T]
h1_1 = 5186*1e-9;                       %[T]
H0 = ( g1_0^2 + g1_1^2 + h1_1^2 )^(1/2);

m = [0.01; 0.05; 0.01]; % parasitic magnetic torque vector: considered a small spacecraft, maximum m (worst case) is [0.1 0.1 0.1]'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Star sensor (mounted on face A4, n4 = [-1; 0; 0])

% Sensor chosen: "Sagitta Star Tracker"
% FOV = 25.4 deg x 25.4 deg
% Limit tracked magnitude: [0:1]
% Accuracy:
    cross_boresight = 2; % ["]
    around_boresight = 10;  % ["]

    cross_boresight = deg2rad(dms2degrees([0,0,cross_boresight])); % [rad]
    around_boresight = deg2rad(dms2degrees([0,0,around_boresight])); % [rad]

    % SOURCE: https://www.cubesatshop.com/product/sagitta-star-tracker/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Catalogue of stars (ref J2000 epoch, sky seen on 2 Jan 2023 from equator above Kourou (0 deg, 0' 0.0''; W 52deg 39' 01'')):

% Comment: take 6 stars spread along different regions of the sky so as to cover the whole orbit
%
% Altair: 19h 50m 46.51s; 8deg 52' 14.3''; Mag: 0.75; 16.73 ly 
%       WARNING: Altair is quite close to the Sun!
% Spica: 13h 25m 11.57s; -11deg 09' 40.75''; Mag: 0.95; 250 ly
% Procyon:  7h 39m 18,1183s; 5deg 13′ 29,975″; Mag: 0.40; 11.4 ly
% Rigel: 05h 14m 32,30  s; -08deg 12′ 06,89″; Mag: 0.15; 860 ly
% Betelgeuse: 05h 55m 10,3053s; 07deg 24′ 25,426″; Mag: 0.45; 500 ly
% Aldebaran: 4h 35m 55,239s; 16deg 30′ 33,488″; Mag: 0.85; 66.64 ly
%       Comment: No bright star between R.A. 20h and 4h along the equator
% Discarded stars (too low magnitude):
% Vega: 18h 36m 56,3364s; 38deg 47′ 1,291″; Mag 0.03; 25.3 ly
% Capella: 5h 16m 41,359s; +45deg 59′ 52,768″; Mag: 0.08; 43 ly
% Rigil Kentaurus: 14h 39m 36,495s; −60deg 50′ 2,308″; Mag: -0.27; 4.364 ly
% Canopus: 06h 23m 57.10988s; 52deg 41′ 44.3810″; Mag: -0.74; 320 ly
% Arcturus: 14h 15m 39.7s; 19deg 10′ 56″; Mag: -0.05; 36.2 ly

% SOURCE: Stellarium
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
% Error matrix is built in Simulink: rotation matrix A312

%       Random number block parameters:
%           Mean: 0
%           Variance: 0.16
%           Seed: 0
%           Sample time: 0.1 (fixed)

% A_epsilon = [1, around_boresight/2 + around_boresight/2*time_effect, -(cross_boresight/2 + cross_boresight/2*time_effect); ...
%              -(around_boresight/2 + around_boresight/2*time_effect), 1, cross_boresight/2 + cross_boresight/2*time_effect; ...
%              cross_boresight/2 + cross_boresight/2*time_effect, -(cross_boresight/2 + cross_boresight/2*time_effect), 1]; % [rad]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Attitude determination

weight = 1/6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pointing

% Nadir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sun sensor (must be installed on same face of star sensor : A4, n4 = [-1; 0; 0])

% I can add a tolerance to the angle at wich the Sun is perceived to turn
% off the star sensor only when it point directly towards the Sun, not just
% when the surface A4 (surface on which the star sensor is mounted) gets exposed to
% the solar radiation.

% Star sensor has a baffle sun rejection angle of 40 deg SOURE: https://www.cubesatshop.com/product/sagitta-star-tracker/

star_sensor_tolerance = 40; % [deg]

star_sensor_tolerance = deg2rad(star_sensor_tolerance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detumbling (magnetic torquer, B-dot method)

Kb = 1e10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control LQR

% Standard state-space form of linear system:
% x_dot = A*x + B*u && y = C*x + D*u

% Extraction of total (main body + deployed elements) inertia moments:
Ix = I(1,1);
Iy = I(2,2);
Iz = I(3,3);

% Construction of matrix A:
A = zeros(length(w_0));
A(1,2) = -((Iz - Iy)/Ix)*n;
A(2,1) = -((Ix - Iz)/Iy)*n;

% Construction of matrix B:
B = diag([1/Ix, 1/Iy, 1/Iz]);

% A, B matrix -> SOURCE: slide 2 of "13b - state space models for attitude dynamics.pdf"

% Construction of matrix Q: weight matrix of the states
Q = zeros(length(w_0));
Q(1,1) = 1/(0.01^2);
Q(2,2) = 1/(0.01^2);
Q(3,3) = 1/(0.001^2);

% Construction of matrix R: weight matrix of the control input
R = zeros(length(w_0));
R(1,1) = 1/((1e-3)^2);
R(2,2) = 1/((1e-3)^2);
R(3,3) = 1/((1e-6)^2);

[K,S,CLP] = lqr(A,B,Q,R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Observer

% C = zeros(3);
% 
% D = zeros(3);
% 
% poles_observer = [-2.25e-06 -0.009 -0.012];
% 
% L_t = place(A', C', poles_observer);
% L = L_t';
% 
% A_hat = A - L*C;
% B_hat = B - L*D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slew Manoeuvre

Kp_slew = [5e-4; 5e-4; 5e-4];
Kd_slew = [5e-5; 5e-5; 5e-5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tracking

Kp_tracking = [5e-2; 5e-2; 5e-2];
Kd_tracking = [5e-3; 0.8*(5e-3); 5e-3];
w_target = [0; 0; n];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Actuators

actuators_matrix = [1 0 0; 0 1 0; 0 0 1];

% Inertia wheel on y (minimum inertial axis), reaction wheel on x and z (respectively maximum inertial axis and intermediate inertial axis)

%% Simulation

% simulation = sim("Model.slx");
% 
% % Attitude
% A_target = simulation.target;
% A_real = simulation.ABN_sensor;
% 
% % External Torques
% torque_srp = simulation.srp;
% torque_gg = simulation.gg;
% torque_mag = simulation.mag;
% 
% % Pointing Error
% pointing_error = simulation.err;
% 
% % Detumbling
% magnetic_torquer_control_torque = simulation.magnetic_torquer;
% 
% % Control Torque
% tracking_torque = simulation.tracking;
% slew_torque = simulation.slew;
% control_torque = simulation.ctrl;
% 
% % RWs and IW torque
% RW_torque = simulation.actuators;
% 
% % Total actuators torque (RW/magnetic torquer)
% total_actuators_torque = simulation.total_actuators_torque; 