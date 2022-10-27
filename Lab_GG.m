
%% Task1
clc
clear all
close all

% Principal moments of inertia: case 1 
Ix = 0.04; % Kgm^2 
Iy = 0.06; % Kgm^2
Iz = 0.08; % Kgm^2
I = [Ix 0 0; 0 Iy 0; 0 0 Iz]; % Inertia matrix
I_inv = inv(I); % Matrice inversa di I

mu_E = astroConstants(13); % Gravitational constant Earth
R_E = astroConstants(23); %Earth's radius
r_orbit = 400; %km
n = sqrt(mu_E/(R_E+ r_orbit)^3); % mean motion
% Initial conditions omega

wx_0 = 0; % Rad/s
wy_0 = 0; % Rad/s
wz_0 = n; % Rad/s


v1 = [0 1 -1]';
v2 = [-1 0 1]';
v3 = [1 -1 0]';

A0 =  eye(3);

w_0 = [wx_0 wy_0 wz_0]'; % vettore condizioni iniziali

w_L = [0 0 n]';

simulation = sim("Lab6_task1");

A_body = simulation.A_body;
omega_b = simulation.w;
omegab_l = simulation.wb_l;

%% caso 2
clc
clear all

% Principal moments of inertia: case 2
Ix = 0.06; % Kgm^2 
Iy = 0.08; % Kgm^2
Iz = 0.04; % Kgm^2
I = [Ix 0 0; 0 Iy 0; 0 0 Iz]; % Inertia matrix
I_inv = inv(I); % Matrice inversa di I

mu_E = astroConstants(13); % Gravitational constant Earth
R_E = astroConstants(23); %Earth's radius
r_orbit = 400; %km
n = sqrt(mu_E/(R_E+ r_orbit)^3); % mean motion
% Initial conditions omega

wx_0 = 1e-06; % Rad/s
wy_0 = 1e-06; % Rad/s
wz_0 = n; % Rad/s


v1 = [0 1 -1]';
v2 = [-1 0 1]';
v3 = [1 -1 0]';

A0 =  eye(3);

w_0 = [wx_0 wy_0 wz_0]'; % vettore condizioni iniziali

w_L = [0 0 n]';

simulation = sim("Lab6_task1");

A_body = simulation.A_body;
omega_b = simulation.w;
omegab_l = simulation.wb_l;

%% Task2 
clc
clear all

% Principal moments of inertia: case 1 
Ix = 0.04; % Kgm^2 
Iy = 0.06; % Kgm^2
Iz = 0.08; % Kgm^2
I = [Ix 0 0; 0 Iy 0; 0 0 Iz]; % Inertia matrix
I_inv = inv(I); % Matrice inversa di I

mu_E = astroConstants(13); % Gravitational constant Earth
R_E = astroConstants(23); %Earth's radius
r_orbit = 400; %km
n = sqrt(mu_E/(R_E+ r_orbit)^3); % mean motion
% Initial conditions omega

wx_0 = 0; % Rad/s
wy_0 = 0; % Rad/s
wz_0 = n; % Rad/s


v1 = [0 1 -1]';
v2 = [-1 0 1]';
v3 = [1 -1 0]';

A0 =  eye(3);

w_0 = [wx_0 wy_0 wz_0]'; % vettore condizioni iniziali

w_L = [0 0 n]';


simulation = sim("Lab6_task2");

A_body = simulation.A_body;
omega_b = simulation.w;
omegab_l = simulation.wb_l;