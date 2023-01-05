clear al
clc
close all

%% Constants:

mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(23);
w_E = deg2rad(15.04/3600);

%% Orbital data:

a = R_E + 550; % Orbital semi-major axis [km]
e = 0.0001824; % Orbital eccentricity
i = deg2rad(53.0545); % Orbital inclination
om = deg2rad(45.9004); % Argument of pericenter
theta0 = deg2rad(0); % True anomaly
n = sqrt(mu_E/(a^3)); % Angular velocity
T = 2*pi/n; % Orbital period

%% coverage:

phi_max = rad2deg(acos(R_E/(a*(1+e))));

%% Orbit propagations:

[r1, v1] = kep2car(a, e, i, deg2rad(10), om, deg2rad(0), mu_E);
tspan = linspace( 0, T, 1000 );

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ ~, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, [r1; v1], options );

figure(1);
earth_sphere('km');
hold on;
plot3( Y1(:,1), Y1(:,2), Y1(:,3), '-r', 'LineWidth', 2);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Orbits');
legend('','Orbit');
axis equal;
grid on;