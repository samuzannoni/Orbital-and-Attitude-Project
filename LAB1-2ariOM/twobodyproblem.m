%% EX1
clear all
close all
clc
% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
% Initial condition
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0, v0];

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T1, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% T è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results

figure()
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;
hold on;
Terra3d

%%
r = [Y(:,1) Y(:,2) Y(:,3)];
v = [Y(:,4) Y(:,5) Y(:,6)];

h = cross(r,v); % angolar momentum
hnorm = vecnorm(h,2,2);

e = 1/mu_E*cross(v,h)-r./vecnorm(r,2,2); % eccentricity
enorm = vecnorm(e,2,2);

rnorm = vecnorm(r,2,2);
vnorm = vecnorm(v,2,2);

ur= r./rnorm;
uh = h./hnorm;
ut = cross(uh,ur);
vr = dot(v,ur,2); % radial velocity
vt = dot(v,ut,2); % transversal velocity
perp= dot(e,h,2); % perpendicolarità

epsilon = vnorm.^2/2-mu_E./rnorm; % specific energy
epsilon_esatta= -mu_E/(2*a); % specific energy 

t_ = linspace(0,5*T,1000);
%% Plot angular momentum
figure(1)
plot(t_,h(:,1));
hold on;
plot(t_,h(:,2));
plot(t_,h(:,3));
plot(t_,hnorm);
grid on;
xlabel('t [s]');
ylabel('h_x, h_y, h_z, ||h|| [km^2/s]');
legend('h_x','h_y','h_z','||h||');
title('angolar momentum');
%% Plot eccentricity
figure(2)
plot(t_,e(:,1));
hold on;
plot(t_,e(:,2));
plot(t_,e(:,3));
plot(t_,enorm);
grid on;
xlabel('t [s]');
ylabel('e_x, e_y, e_z, ||e|| [-]');
legend('e_x','e_y','e_z','||e||');
title('eccentricity vector');
%% Plot radial and transversal velocity
figure(3)
plot(t_,vr);
hold on;
plot(t_,vt);
grid on;
xlabel('t [s]');
ylabel('v_x, v_y [km/s]');
title('radial and transversal velocity');
legend('v_r','v_t');
%% Plot e-h dot
figure(4)
plot(t_,perp);
grid on;
xlabel('t [s]');
ylabel('|e * h| [ km^2/s]');
title('e-h dot product');

%% Plot specific energy
figure(5)
plot(t_,epsilon);
grid on;
xlabel('t [s]');
ylabel('epsilon [ km^2/s^2]');
title('energy');



%% HIGHLY ECCENTRIC AND INCLINED ORBIT
clear all
close all
clc
% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]

% Initial Conditions
r0 = [ 6495; -970; -3622 ]; %[km]
v0 = [ 4.752; 2.130; 7.950 ]; %[km/s]
y0 = [r0, v0];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T2, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% T è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results
figure()
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;
hold on;
Terra3d;
%%
r = [Y(:,1) Y(:,2) Y(:,3)];
v = [Y(:,4) Y(:,5) Y(:,6)];

h = cross(r,v); % angolar momentum
hnorm = vecnorm(h,2,2);

rnorm = vecnorm(r,2,2);
vnorm= vecnorm(v,2,2);

e = 1/mu_E*cross(v,h)-r./rnorm; % eccentricity
enorm = vecnorm(e,2,2);

ur= r./rnorm;
uh = h./hnorm;
ut = cross(uh,ur);
vr = dot(v,ur,2); % radial velocity
vt = dot(v,ut,2); % transversal velocity
perp= dot(e,h,2); % perpendicolarità

epsilon = vnorm.^2/2-mu_E./rnorm; % specific energy
epsilon_esatta= -mu_E/(2*a); % specific energy 

t_ = linspace(0,5*T,1000);
%% Plot angular momentum
figure(1)
plot(t_,h(:,1));
hold on;
plot(t_,h(:,2));
plot(t_,h(:,3));
plot(t_,hnorm);
grid on;
xlabel('t [s]');
ylabel('h_x, h_y, h_z, ||h|| [km^2/s]');
legend('h_x','h_y','h_z','||h||');
title('angolar momentum');
%% Plot eccentricity
figure(2)
plot(t_,e(:,1));
hold on;
plot(t_,e(:,2));
plot(t_,e(:,3));
plot(t_,enorm);
grid on;
xlabel('t [s]');
ylabel('e_x, e_y, e_z, ||e|| [-]');
legend('e_x','e_y','e_z','||e||');
title('eccentricity vector');
%% Plot radial and transversal velocity
figure(3)
plot(t_,vr);
hold on;
plot(t_,vt);
grid on;
xlabel('t [s]');
ylabel('v_x, v_y [km/s]');
title('radial and transversal velocity');
legend('v_r','v_t');
%% Plot e-h dot
figure(4)
plot(t_,perp);
grid on;
xlabel('t [s]');
ylabel('|e * h| [ km^2/s]');
title('e-h dot product');

%% Plot specific energy
figure(5)
plot(t_,epsilon);
grid on;
xlabel('t [s]');
ylabel('epsilon [ km^2/s^2]');
title('energy');


