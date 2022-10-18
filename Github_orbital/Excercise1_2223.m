%%
clear all
close all
clc
% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
% Initial condition
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T1, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% T1 è tspan in colonna
% Y sono le componenti posizioni e velociutà
% Plot the results
figure(1)
Terra3d
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on

%%
r = [Y(:,1) Y(:,2) Y(:,3)];
v = [Y(:,4) Y(:,5) Y(:,6)];
h = cross(r,v);

e = 1/mu_E*cross(v,h)-r./vecnorm(r,2,2);
ur= r./vecnorm(r,2,2);
uh = h./vecnorm(h,2,2);
ut = cross(uh,ur);
vr = dot(v',ur');
vt = dot(v',ut');
perp= dot(e',h'); % perpendicolarità tra e e h

epsilon = (vecnorm(v,2,2)).^2/2-mu_E./vecnorm(r,2,2)
epsilon_esatta= -mu_E/(2*a); % energia specifica

figure(2)
t_ = linspace(0,5*T,1000);
plot(t_,h(:,1));
hold on;
plot(t_,h(:,2));
plot(t_,h(:,3));
grid on;
title('Momento angolare')
xlabel('t [s]')
ylabel('h [km^2/s]')
legend('h_x','h_y','h_z')

figure(3)
plot(t_,e(:,1));
hold on;
plot(t_,e(:,2));
plot(t_,e(:,3));
title('eccentricity')
xlabel('t [s]')
ylabel('e')
legend('e_x','e_y','e_z')
grid on

figure(4)
plot(t_,vr);
hold on;
plot(t_,vt);
grid on;
legend('v_r','v_t')
xlabel('t [s]')
ylabel('v [km/s]')
title(' Radial and transversal velocity')

figure(5)
plot(t_,perp);
grid on;
title('e-h product')


figure(6)
plot(t_,epsilon);
grid on;
title('Energy')

%% Highly eccentric and inclined orbit
clear all
close all
clc
% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
% Initial condition
r0 = [ 6495; -970; 0-3622]; % [km]
v0 = [ 4.752; 2.130; 7.950 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T1, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% T1 è tspan in colonna
% Y sono le componenti posizioni e velociutà
% Plot the results
figure(1)
Terra3d
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on

r = [Y(:,1) Y(:,2) Y(:,3)];
v = [Y(:,4) Y(:,5) Y(:,6)];
h = cross(r,v);

e = 1/mu_E*cross(v,h)-r./vecnorm(r,2,2);
ur= r./vecnorm(r,2,2);
uh = h./vecnorm(h,2,2);
ut = cross(uh,ur);
vr = dot(v',ur');
vt = dot(v',ut');
perp= dot(e',h'); % perpendicolarità tra e e h

epsilon = (vecnorm(v,2,2)).^2/2-mu_E./vecnorm(r,2,2)
epsilon_esatta= -mu_E/(2*a); % energia specifica

figure(2)
t_ = linspace(0,5*T,1000);
plot(t_,h(:,1));
hold on;
plot(t_,h(:,2));
plot(t_,h(:,3));
grid on;
title('Momento angolare')
xlabel('t [s]')
ylabel('h [km^2/s]')
legend('h_x','h_y','h_z')

figure(3)
plot(t_,e(:,1));
hold on;
plot(t_,e(:,2));
plot(t_,e(:,3));
title('eccentricity')
xlabel('t [s]')
ylabel('e')
legend('e_x','e_y','e_z')
grid on

figure(4)
plot(t_,vr);
hold on;
plot(t_,vt);
grid on;
legend('v_r','v_t')
xlabel('t [s]')
ylabel('v [km/s]')
title(' Radial and transversal velocity')

figure(5)
plot(t_,perp);
grid on;
title('e-h product')


figure(6)
plot(t_,epsilon);
grid on;
title('Energy') %commento
