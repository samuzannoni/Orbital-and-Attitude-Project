clear al
clc
close all

%% Constants:

mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(23);
w_E = deg2rad(15.04/3600);

%% Orbital data:

i = deg2rad(50);
T = 10.8*3600;
e = 0.7;
om = deg2rad(270);
a = ((T/(2*pi))^2*mu_E)^(1/3);

%% Orbit propagations:

[r1, v1] = kep2car(a, e, i, deg2rad(10), om, deg2rad(0), mu_E);

tspan = linspace( 0, 3*T, 1000 );

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ ~, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, [r1; v1], options );

%% Plot Orbit:

figure(1);
earth_sphere('km');
hold on;
plot3( Y1(:,1), Y1(:,2), Y1(:,3), '-y');

xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Orbits');
legend('','Orbit1', 'Orbit2', 'Orbit3');
axis equal;
grid on;

%% Compute Groundtrack:

theta_G0_1 = zeros(1, length(tspan)); % Greenwhich longitude
t0_1 = zeros(1, length(tspan));

[alpha1, delta1, lon1, lat1]=groundTrack(r1, v1, theta_G0_1, tspan, w_E, mu_E, t0_1);


lat1=rad2deg(lat1);
lon1=wrapTo180(rad2deg(lon1));

% Control to avoid discontinuities:
for z = 1:length(lon1)-1
    if (lon1(z)>0 && lon1(z+1)<0)
        lon1(z) = NaN;
        lat1(z) = NaN;
    end
end

%% Plot Groundtrack:

figure(2);
planisphere = imread("earth_planisphere.jpg");
image('CData', planisphere, 'XData', [-180 180], 'YData', [90 -90]);
hold on
plot(lon1, lat1, '-y');
axis([-180 180 -90 90]);
plot(lon1(1),lat1(1),'oy',lon1(end),lat1(end),'sy','MarkerSize',10, 'LineWidth',2); % insert marker

title('Groundtrack');
legend('Groundtrack1', 'Start', 'End');
xlabel('Longitude');
ylabel('Latitude');
