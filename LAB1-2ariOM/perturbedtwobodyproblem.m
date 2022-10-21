%% EX2 Quasi_CIRCULAR ORBIT 5T-PERIOD
clear all
close all
clc
% Physical parameters
out = astroConstants([13,23,9]); % Earth's gravitational parameter [km^3/s^2]
% Initial condition
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0, v0];
mu_E = out(1);
Re = out(2);
J2 = out(3);
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ Tp, Yp] = ode113( @(t,y) ode_2bpp(t,y,mu_E, Re, J2), tspan, y0, options );
[ T1, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% T è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results

figure()
plot3( Yp(:,1), Yp(:,2), Yp(:,3), '-' );
%hold on;
%plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Perturbed-two-body problem orbit');
%legend( 'perturbed orbit','unperturbed orbit');
axis equal;
grid on;
hold on;
Terra3d

%% Calulus
% Calculate parameters of perturbed 2bp
rp = [Yp(:,1) Yp(:,2) Yp(:,3)];
vp = [Yp(:,4) Yp(:,5) Yp(:,6)];

hp = cross(rp,vp); % angolar momentum
hnormp= vecnorm(hp,2,2);
rnormp = vecnorm(rp,2,2);
vnormp = vecnorm(vp,2,2);

ep = 1/mu_E*cross(vp,hp)-rp./rnormp; % eccentricity
enormp = vecnorm(ep,2,2);

urp= rp./rnormp;
uhp = hp./hnormp;
utp = cross(uhp,urp);
vrp = dot(vp,urp,2); % radial velocity
vtp = dot(vp,utp,2); % transversal velocity
perpp= dot(ep,hp,2); % perpendicolarità

epsilonp = vnormp.^2/2-mu_E./rnormp; % specific energy
epsilon_esattap= -mu_E/(2*a); % specific energy 


% Calculate parameters of 2bp
r = [Y(:,1) Y(:,2) Y(:,3)];
v = [Y(:,4) Y(:,5) Y(:,6)];

h = cross(r,v); % angolar momentum
hnorm= vecnorm(h,2,2);
rnorm = vecnorm(r,2,2);
vnorm = vecnorm(v,2,2);

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
plot(t_,hp(:,1),'b-.',t_,hp(:,2),'g-.',t_,hp(:,3),'r-.',t_,hnormp,'k-.');
hold on;
plot(t_,h(:,1),'b--',t_,h(:,2),'g--',t_,h(:,3),'r--',t_,hnorm,'k--');
grid on;
xlabel('t [s]');
ylabel('h_x, h_y, h_z, ||h|| [km^2/s]');
legend('h_xJ2perturbed','h_yJ2perturbed','h_zJ2perturbed','||h||J2perturbed','h_x','h_y','h_z','||h||');
title('angolar momentum');
%% Plot eccentricity
figure(2)
plot(t_,ep(:,1),'b-.',t_,ep(:,2),'g-.',t_,ep(:,3),'r-.',t_,enormp,'k-.');
hold on;
plot(t_,e(:,1),'b--',t_,e(:,2),'g--',t_,e(:,3),'r--',t_,enorm,'k--');
grid on;
xlabel('t [s]');
ylabel('e_x, e_y, e_z, ||e||');
legend('e_xJ2perturbed','e_yJ2perturbed','e_zJ2perturbed','||e||J2perturbed','e_x','e_y','e_z','||e||');
title('eccentricity vector');
%% Plot radial and transversal velocity
figure(3)
plot(t_,vrp,'b-.',t_,vtp,'g-.');
hold on;
plot(t_,vr,'b-.',t_,vt,'g-.');
grid on;
title('radial and transversal velocity');
legend('v_rJ2perturbed','v_tJ2perturbed','v_r','v_t');
%% Plot e-h dot
figure(4)
plot(t_,perpp,'r');
hold on;
plot(t_,perp,'b')
grid on;
xlabel('t [s]');
ylabel('|e * h| [ km^2/s]');
title('e dot h');
legend('e-h dot product J2perturbed','e-h dot product');

%% Plot specific energy
figure(5)
plot(t_,epsilonp,'r');
hold on;
plot(t_,epsilon,'b');
grid on;
legend('energy J2perturbed','energy');

%% Year
clear all
close all
clc

% Physical parameters
out = astroConstants([13,23,9]); % Earth's gravitational parameter [km^3/s^2]
% Initial condition
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0, v0];
mu_E = out(1);
Re = out(2);
J2 = out(3);
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
year = 31536000;
tspan = linspace( 0, year, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ Tp, Yp] = ode113( @(t,y) ode_2bpp(t,y,mu_E, Re, J2), tspan, y0, options );
[ T1, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% T è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results

figure()
plot3( Yp(:,1), Yp(:,2), Yp(:,3), '-' );
%hold on;
%plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Perturbed-two-body problem orbit');
%legend( 'perturbed orbit','unperturbed orbit');
axis equal;
grid on;
hold on;
Terra3d

%% Calulus
% Calculate parameters of perturbed 2bp
rp = [Yp(:,1) Yp(:,2) Yp(:,3)];
vp = [Yp(:,4) Yp(:,5) Yp(:,6)];

hp = cross(rp,vp); % angolar momentum
hnormp= vecnorm(hp,2,2);
rnormp = vecnorm(rp,2,2);
vnormp = vecnorm(vp,2,2);

ep = 1/mu_E*cross(vp,hp)-rp./rnormp; % eccentricity
enormp = vecnorm(ep,2,2);

urp= rp./rnormp;
uhp = hp./hnormp;
utp = cross(uhp,urp);
vrp = dot(vp,urp,2); % radial velocity
vtp = dot(vp,utp,2); % transversal velocity
perpp= dot(ep,hp,2); % perpendicolarità

epsilonp = vnormp.^2/2-mu_E./rnormp; % specific energy
epsilon_esattap= -mu_E/(2*a); % specific energy 


% Calculate parameters of 2bp
r = [Y(:,1) Y(:,2) Y(:,3)];
v = [Y(:,4) Y(:,5) Y(:,6)];

h = cross(r,v); % angolar momentum
hnorm= vecnorm(h,2,2);
rnorm = vecnorm(r,2,2);
vnorm = vecnorm(v,2,2);

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

t_ = linspace(0,year,1000);
%% Plot angular momentum
figure(1)
plot(t_,hp(:,1),'b-.',t_,hp(:,2),'g-.',t_,hp(:,3),'r-.',t_,hnormp,'k-.');
hold on;
plot(t_,h(:,1),'b--',t_,h(:,2),'g--',t_,h(:,3),'r--',t_,hnorm,'k--');
grid on;
xlabel('t [s]');
ylabel('h_x, h_y, h_z, ||h|| [km^2/s]');
legend('h_xJ2perturbed','h_yJ2perturbed','h_zJ2perturbed','||h||J2perturbed','h_x','h_y','h_z','||h||');
title('angolar momentum');
%% Plot eccentricity
figure(2)
plot(t_,ep(:,1),'b-.',t_,ep(:,2),'g-.',t_,ep(:,3),'r-.',t_,enormp,'k-.');
hold on;
plot(t_,e(:,1),'b--',t_,e(:,2),'g--',t_,e(:,3),'r--',t_,enorm,'k--');
grid on;
xlabel('t [s]');
ylabel('e_x, e_y, e_z, ||e||');
legend('e_xJ2perturbed','e_yJ2perturbed','e_zJ2perturbed','||e||J2perturbed','e_x','e_y','e_z','||e||');
title('eccentricity vector');
%% Plot radial and transversal velocity
figure(3)
plot(t_,vrp,'b-.',t_,vtp,'g-.');
hold on;
plot(t_,vr,'b-.',t_,vt,'g-.');
grid on;
title('radial and transversal velocity');
legend('v_rJ2perturbed','v_tJ2perturbed','v_r','v_t');
%% Plot e-h dot
figure(4)
plot(t_,perpp,'r');
hold on;
plot(t_,perp,'b')
grid on;
xlabel('t [s]');
ylabel('|e * h| [ km^2/s]');
title('e dot h');
legend('e-h dot product J2perturbed','e-h dot product');

%% Plot specific energy
figure(5)
plot(t_,epsilonp,'r');
hold on;
plot(t_,epsilon,'b');
grid on;
legend('energy J2perturbed','energy');


%% HIGHLY ECCENTRIC AND INCLINED ORBIT 5T-PERIOD
clear all
close all
clc
% Physical parameters
out = astroConstants([13,23,9]); % Earth's gravitational parameter [km^3/s^2]
% Initial condition
r0 = [ 6495; -970; -3622 ]; %[km]
v0 = [ 4.752; 2.130; 7.950 ]; %[km/s]
y0 = [r0, v0];
mu_E = out(1);
Re = out(2);
J2 = out(3);
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ Tp, Yp] = ode113( @(t,y) ode_2bpp(t,y,mu_E, Re, J2), tspan, y0, options );
[ T1, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% T è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results

figure()
plot3( Yp(:,1), Yp(:,2), Yp(:,3), '-' );
%hold on;
%plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Perturbed-two-body problem orbit');
%legend( 'perturbed orbit','unperturbed orbit');
axis equal;
grid on;
hold on;
Terra3d

%% Calulus
% Calculate parameters of perturbed 2bp
rp = [Yp(:,1) Yp(:,2) Yp(:,3)];
vp = [Yp(:,4) Yp(:,5) Yp(:,6)];

hp = cross(rp,vp); % angolar momentum
hnormp= vecnorm(hp,2,2);
rnormp = vecnorm(rp,2,2);
vnormp = vecnorm(vp,2,2);

ep = 1/mu_E*cross(vp,hp)-rp./rnormp; % eccentricity
enormp = vecnorm(ep,2,2);

urp= rp./rnormp;
uhp = hp./hnormp;
utp = cross(uhp,urp);
vrp = dot(vp,urp,2); % radial velocity
vtp = dot(vp,utp,2); % transversal velocity
perpp= dot(ep,hp,2); % perpendicolarità

epsilonp = vnormp.^2/2-mu_E./rnormp; % specific energy
epsilon_esattap= -mu_E/(2*a); % specific energy 


% Calculate parameters of 2bp
r = [Y(:,1) Y(:,2) Y(:,3)];
v = [Y(:,4) Y(:,5) Y(:,6)];

h = cross(r,v); % angolar momentum
hnorm= vecnorm(h,2,2);
rnorm = vecnorm(r,2,2);
vnorm = vecnorm(v,2,2);

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
plot(t_,hp(:,1),'b-.',t_,hp(:,2),'g-.',t_,hp(:,3),'r-.',t_,hnormp,'k-.');
hold on;
plot(t_,h(:,1),'b--',t_,h(:,2),'g--',t_,h(:,3),'r--',t_,hnorm,'k--');
grid on;
xlabel('t [s]');
ylabel('h_x, h_y, h_z, ||h|| [km^2/s]');
legend('h_xJ2perturbed','h_yJ2perturbed','h_zJ2perturbed','||h||J2perturbed','h_x','h_y','h_z','||h||');
title('angolar momentum');
%% Plot eccentricity
figure(2)
plot(t_,ep(:,1),'b-.',t_,ep(:,2),'g-.',t_,ep(:,3),'r-.',t_,enormp,'k-.');
hold on;
plot(t_,e(:,1),'b--',t_,e(:,2),'g--',t_,e(:,3),'r--',t_,enorm,'k--');
grid on;
xlabel('t [s]');
ylabel('e_x, e_y, e_z, ||e||');
legend('e_xJ2perturbed','e_yJ2perturbed','e_zJ2perturbed','||e||J2perturbed','e_x','e_y','e_z','||e||');
title('eccentricity vector');
%% Plot radial and transversal velocity
figure(3)
plot(t_,vrp,'b-.',t_,vtp,'g-.');
hold on;
plot(t_,vr,'b-.',t_,vt,'g-.');
grid on;
title('radial and transversal velocity');
legend('v_rJ2perturbed','v_tJ2perturbed','v_r','v_t');
%% Plot e-h dot
figure(4)
plot(t_,perpp,'r');
hold on;
plot(t_,perp,'b')
grid on;
xlabel('t [s]');
ylabel('|e * h| [ km^2/s]');
title('e dot h');
legend('e-h dot product J2perturbed','e-h dot product');

%% Plot specific energy
figure(5)
plot(t_,epsilonp,'r');
hold on;
plot(t_,epsilon,'b');
grid on;
legend('energy J2perturbed','energy');



%% YEAR
clear all
close all
clc
% Physical parameters
out = astroConstants([13,23,9]); % Earth's gravitational parameter [km^3/s^2]
% Initial condition
r0 = [ 6495; -970; -3622 ]; %[km]
v0 = [ 4.752; 2.130; 7.950 ]; %[km/s]
y0 = [r0, v0];
mu_E = out(1);
Re = out(2);
J2 = out(3);
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
year = 31536000;
i = year/T;
tspan = linspace( 0, i*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ Tp, Yp] = ode113( @(t,y) ode_2bpp(t,y,mu_E, Re, J2), tspan, y0, options );
[ T1, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
% T è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results

figure()
plot3( Yp(:,1), Yp(:,2), Yp(:,3), '-' );
%hold on;
%plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Perturbed-two-body problem orbit');
%legend( 'perturbed orbit','unperturbed orbit');
axis equal;
grid on;
hold on;
Terra3d

%% Calulus
% Calculate parameters of perturbed 2bp
rp = [Yp(:,1) Yp(:,2) Yp(:,3)];
vp = [Yp(:,4) Yp(:,5) Yp(:,6)];

hp = cross(rp,vp); % angolar momentum
hnormp= vecnorm(hp,2,2);
rnormp = vecnorm(rp,2,2);
vnormp = vecnorm(vp,2,2);

ep = 1/mu_E*cross(vp,hp)-rp./rnormp; % eccentricity
enormp = vecnorm(ep,2,2);

urp= rp./rnormp;
uhp = hp./hnormp;
utp = cross(uhp,urp);
vrp = dot(vp,urp,2); % radial velocity
vtp = dot(vp,utp,2); % transversal velocity
perpp= dot(ep,hp,2); % perpendicolarità

epsilonp = vnormp.^2/2-mu_E./rnormp; % specific energy
epsilon_esattap= -mu_E/(2*a); % specific energy 


% Calculate parameters of 2bp
r = [Y(:,1) Y(:,2) Y(:,3)];
v = [Y(:,4) Y(:,5) Y(:,6)];

h = cross(r,v); % angolar momentum
hnorm= vecnorm(h,2,2);
rnorm = vecnorm(r,2,2);
vnorm = vecnorm(v,2,2);

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

t_ = linspace(0,i*T,1000);
%% Plot angular momentum
figure(1)
plot(t_,hp(:,1),'b-.',t_,hp(:,2),'g-.',t_,hp(:,3),'r-.',t_,hnormp,'k-.');
hold on;
plot(t_,h(:,1),'b--',t_,h(:,2),'g--',t_,h(:,3),'r--',t_,hnorm,'k--');
grid on;
xlabel('t [s]');
ylabel('h_x, h_y, h_z, ||h|| [km^2/s]');
legend('h_xJ2perturbed','h_yJ2perturbed','h_zJ2perturbed','||h||J2perturbed','h_x','h_y','h_z','||h||');
title('angolar momentum');
%% Plot eccentricity
figure(2)
plot(t_,ep(:,1),'b-.',t_,ep(:,2),'g-.',t_,ep(:,3),'r-.',t_,enormp,'k-.');
hold on;
plot(t_,e(:,1),'b--',t_,e(:,2),'g--',t_,e(:,3),'r--',t_,enorm,'k--');
grid on;
xlabel('t [s]');
ylabel('e_x, e_y, e_z, ||e||');
legend('e_xJ2perturbed','e_yJ2perturbed','e_zJ2perturbed','||e||J2perturbed','e_x','e_y','e_z','||e||');
title('eccentricity vector');
%% Plot radial and transversal velocity
figure(3)
plot(t_,vrp,'b-.',t_,vtp,'g-.');
hold on;
plot(t_,vr,'b-.',t_,vt,'g-.');
grid on;
title('radial and transversal velocity');
legend('v_rJ2perturbed','v_tJ2perturbed','v_r','v_t');
%% Plot e-h dot
figure(4)
plot(t_,perpp,'r');
hold on;
plot(t_,perp,'b')
grid on;
xlabel('t [s]');
ylabel('|e * h| [ km^2/s]');
title('e dot h');
legend('e-h dot product J2perturbed','e-h dot product');

%% Plot specific energy
figure(5)
plot(t_,epsilonp,'r');
hold on;
plot(t_,epsilon,'b');
grid on;
legend('energy J2perturbed','energy');


