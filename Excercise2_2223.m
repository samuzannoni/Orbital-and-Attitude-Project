clc
clear all 
close all


% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_Earth = astroConstants(23); % Radius of the Earth
J2 = astroConstants(9); % Coefficient of oblateness perturbation


% Initial condition, orbit not inclined
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
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem

% T1, T2 è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results over 5 periods
figure(1)
Terra3d
plot3( Y1(:,1), Y1(:,2), Y1(:,3), '--k')

view([-40 20])
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
hold on
plot3( Y2(:,1), Y2(:,2), Y2(:,3), '--r' )
legend('Unperturbed orbit','Perturbed orbit')
grid on


%% Comparison unperturbed orbit vs perturbed orbit over 1 year period

tspan = linspace( 0, 24*3600*365, 100000 ); %increased the point of tspan
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem

% T1, T2 è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results over 5 periods
figure(2)
%N_orbit = T2./T;
Terra3d
cb = colorbar();
plot3( Y1(:,1), Y1(:,2), Y1(:,3), '--b')
view([-45 20])
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
hold on
plot3( Y2(:,1), Y2(:,2), Y2(:,3), '--y' )
legend('Unperturbed orbit','Perturbed orbit (J2)')
grid on

%% angular momentum
clc
close all

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


% Comparison between angular momentum of unperturbed vs perturbed orbit
% over 5 periods

r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);

figure(3)

plot(tspan,h1(:,1),'--b','LineWidth',1.5);
hold on;
plot(tspan,h1(:,2),'--r','LineWidth',1.5);
plot(tspan,h1(:,3),'--g','LineWidth',1.5);
hold on
plot(tspan,h2_norm,'--k','Linewidth',2)
grid on;
title('Momento angolare 5 periods')
xlabel('t [s]')
ylabel('h [km^2/s]')


hold on
plot(tspan,h2(:,1),'-.b','LineWidth',1.5);
hold on;
plot(tspan,h2(:,2),'-.r','LineWidth',1.5);
plot(tspan,h2(:,3),'-.g','LineWidth',1.5);
hold on
plot(tspan,h2_norm,'-.k','Linewidth',1.5)

legend('h_x, 2BP','h_y,2BP','h_z,2BP','h_x, J2 perturbed','h_y,J2 Perturbed','h_z,J2 Perturbed')


%% Comparison between angular momentum of unperturbed vs perturbed orbit
%% over 1 year
clc
close all

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 24*3600*365, 10000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);

figure(4)

plot(tspan,h1(:,1),'--b','LineWidth',1.5);
hold on;
plot(tspan,h1(:,2),'--r','LineWidth',1.5);
plot(tspan,h1(:,3),'--g','LineWidth',1.5);
hold on
plot(tspan,h1_norm,'--k','Linewidth',2)
grid on;
title('Momento angolare 1 year')
xlabel('t [s]')
ylabel('h [km^2/s]')


hold on
plot(tspan,h2(:,1),'-.b','LineWidth',1.5);
hold on;
plot(tspan,h2(:,2),'-.r','LineWidth',1.5);
plot(tspan,h2(:,3),'-.g','LineWidth',1.5);
hold on
plot(tspan,h2_norm,'-.k','Linewidth',1.5)

legend('h_x, 2BP','h_y,2BP','h_z,2BP','h_x, J2 perturbed','h_y,J2 Perturbed','h_z,J2 Perturbed')

%% Eccentricity vector 5 periods
clc
clear all 
close all

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_Earth = astroConstants(23); % Radius of the Earth
J2 = astroConstants(9); % Coefficient of oblateness perturbation


% Initial condition, orbit not inclined
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 )';
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);

e1 = 1/mu_E*cross(v1,h1)-r1./vecnorm(r1,2,2); % eccentricity vector unperturbed orbit
e1_norm = vecnorm(e1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);

e2 = 1/mu_E*cross(v2,h2)-r2./vecnorm(r2,2,2); % eccentricity vector perturbed orbit
e2_norm = vecnorm(e2,2,2);

figure(4)
plot(tspan,e1(:,1),'--b',tspan,e1(:,2),'--r',tspan,e1(:,3),'--g',tspan,e1_norm,'--k','LineWidth',1.5)
grid on

hold on
plot(tspan,e2(:,1),'--b',tspan,e2(:,2),'--r',tspan,e2(:,3),'--g',tspan,e2_norm,'--k','LineWidth',1.5)

title('Eccentricity vector - 5 periods')
legend('e_x, 2BP','e_y, 2BP','e_z, 2BP','|e|','e_x, J2 perturbed','e_y, J2 perturbed','e_z, J2 perturbed','|e|, J2 perturbed')

%% Eccentricity vector 1 year period
clc
clear all 
close all

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_Earth = astroConstants(23); % Radius of the Earth
J2 = astroConstants(9); % Coefficient of oblateness perturbation


% Initial condition, orbit not inclined
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 24*3600*365, 10000 )';
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);

e1 = 1/mu_E*cross(v1,h1)-r1./vecnorm(r1,2,2); % eccentricity vector unperturbed orbit
e1_norm = vecnorm(e1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);

e2 = 1/mu_E*cross(v2,h2)-r2./vecnorm(r2,2,2); % eccentricity vector perturbed orbit
e2_norm = vecnorm(e2,2,2);

figure(4)
plot(tspan,e1(:,1),'--b',tspan,e1(:,2),'--r',tspan,e1(:,3),'--g',tspan,e1_norm,'--k','LineWidth',1.5)
grid on

hold on
plot(tspan,e2(:,1),'--b',tspan,e2(:,2),'--r',tspan,e2(:,3),'--g',tspan,e2_norm,'--k','LineWidth',1.5)

title('Eccentricity vector - 5 periods')
legend('e_x, 2BP','e_y, 2BP','e_z, 2BP','|e|','e_x, J2 perturbed','e_y, J2 perturbed','e_z, J2 perturbed','|e|, J2 perturbed')


%% Radial and transversal velocity - 5 periods

clc
clear all 
close all

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_Earth = astroConstants(23); % Radius of the Earth
J2 = astroConstants(9); % Coefficient of oblateness perturbation


% Initial condition, orbit not inclined
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 )';
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);

ur1= r1./vecnorm(r1,2,2);
uh1 = h1./vecnorm(h1,2,2);
ut1 = cross(uh1,ur1);
vr1 = dot(v1',ur1');
vt1 = dot(v1',ut1');

ur2= r2./vecnorm(r2,2,2);
uh2 = h2./vecnorm(h2,2,2);
ut2 = cross(uh2,ur2);
vr2 = dot(v2',ur2');
vt2 = dot(v2',ut2');

figure(5)
plot(tspan,vr1,'--b',tspan,vt1,'--r','LineWidth',1.5)
grid on
hold on
plot(tspan,vr2,':g',tspan,vt2,':k','LineWidth',1.5)
title('Radial and transversal velocity - 5 periods')
xlabel('t [s]')
ylabel('v_r, v_t [km/s]')
legend('v_r, 2BP','v_t, 2BP','v_r, J2 Pertutbed','v_t, J2 perturbed')

%% %% Radial and transversal velocity - 1 year periods

T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 24*3600*365, 10000 )';
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);

ur1= r1./vecnorm(r1,2,2);
uh1 = h1./vecnorm(h1,2,2);
ut1 = cross(uh1,ur1);
vr1 = dot(v1',ur1');
vt1 = dot(v1',ut1');

ur2= r2./vecnorm(r2,2,2);
uh2 = h2./vecnorm(h2,2,2);
ut2 = cross(uh2,ur2);
vr2 = dot(v2',ur2');
vt2 = dot(v2',ut2');

figure(5)
plot(tspan,vr1,'--b',tspan,vt1,'--r','LineWidth',1.5)
grid on
hold on
plot(tspan,vr2,':g',tspan,vt2,':k','LineWidth',1.5)
title('Radial and transversal velocity - 1 year period')
xlabel('t [s]')
ylabel('v_r, v_t [km/s]')
legend('v_r, 2BP','v_t, 2BP','v_r, J2 Pertutbed','v_t, J2 perturbed')

%% Energy over 5 periods
clc
clear all
close all

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_Earth = astroConstants(23); % Radius of the Earth
J2 = astroConstants(9); % Coefficient of oblateness perturbation


% Initial condition, orbit not inclined
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 )';
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);

epsilon1 = (vecnorm(v1,2,2)).^2/2-mu_E./vecnorm(r1,2,2); %2BP
epsilon2 = (vecnorm(v2,2,2)).^2/2-mu_E./vecnorm(r2,2,2); %Perturbed

plot(tspan,epsilon1,'-b',tspan,epsilon2,'--r','LineWidth',1.5)
grid on
xlabel('t [s]')
ylabel('\epsilon [km^2/s^2]')
title('Energy over 5 periods')
legend('2BP','J2 Perturbed')

%% Energy over 1 year period
clc
clear all
close all

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_Earth = astroConstants(23); % Radius of the Earth
J2 = astroConstants(9); % Coefficient of oblateness perturbation


% Initial condition, orbit not inclined
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 24*3600*365, 10000 )';
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);

epsilon1 = (vecnorm(v1,2,2)).^2/2-mu_E./vecnorm(r1,2,2); %2BP
epsilon2 = (vecnorm(v2,2,2)).^2/2-mu_E./vecnorm(r2,2,2); %Perturbed

plot(tspan,epsilon1,'-b',tspan,epsilon2,'--r','LineWidth',1.5)
grid on
xlabel('t [s]')
ylabel('\epsilon [km^2/s^2]')
title('Energy over 1 year period')
legend('2BP','J2 Perturbed')

%% e-h dot product 5 period

clc
clear all
close all

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_Earth = astroConstants(23); % Radius of the Earth
J2 = astroConstants(9); % Coefficient of oblateness perturbation


% Initial condition, orbit not inclined
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 )';
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);


e1 = 1/mu_E*cross(v1,h1)-r1./vecnorm(r1,2,2); % eccentricity vector unperturbed orbit
e1_norm = vecnorm(e1,2,2);

e2 = 1/mu_E*cross(v2,h2)-r2./vecnorm(r2,2,2); % eccentricity vector perturbed orbit
e2_norm = vecnorm(e2,2,2);

perp1= dot(e1',h1'); % perpendicolarità tra e e h
perp2= dot(e2',h2'); % perpendicolarità tra e e h

plot(tspan,perp1,'-b',tspan,perp2,'--r','LineWidth',1.5)
grid on
xlabel('t [s]')
ylabel('|e h| [km^2/s]')
legend('2BP','J2-Perturbed')
title('e-h dot product (5 periods)')

%% %% e-h dot product

clc
clear all
close all

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_Earth = astroConstants(23); % Radius of the Earth
J2 = astroConstants(9); % Coefficient of oblateness perturbation


% Initial condition, orbit not inclined
r0 = [ 26578.137; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 24*3600*365, 10000 )';
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);


e1 = 1/mu_E*cross(v1,h1)-r1./vecnorm(r1,2,2); % eccentricity vector unperturbed orbit
e1_norm = vecnorm(e1,2,2);

e2 = 1/mu_E*cross(v2,h2)-r2./vecnorm(r2,2,2); % eccentricity vector perturbed orbit
e2_norm = vecnorm(e2,2,2);

perp1= dot(e1',h1'); % perpendicolarità tra e e h
perp2= dot(e2',h2'); % perpendicolarità tra e e h

plot(tspan,perp1,'-b',tspan,perp2,'--r','LineWidth',1.5)
grid on
xlabel('t [s]')
ylabel('|e h| [km^2/s]')
legend('2BP','J2-Perturbed')
title('e-h dot product (1 year period)')


%% Perturbed 2 Body problem : Highly eccentric and inclined orbit
clc
close all
% Initial condition, inclined orbit
r0 = [ 6495; -970; -3622 ]; % [km]
v0 = [ 4.752; 2.130; 7.950 ]; % [km/s]
y0 = [ r0; v0 ];
% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem

% T1, T2 è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results over 5 periods
figure(1)
Terra3d
plot3( Y1(:,1), Y1(:,2), Y1(:,3), '--k')

view([-30 20])
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit - 5 orbits');
axis equal;
hold on
plot3( Y2(:,1), Y2(:,2), Y2(:,3), '--r' )
legend('Unperturbed orbit','Perturbed orbit')
grid on

%% Higlhy eccentric and inclined orbit - 1 year period
tspan = linspace( 0, 3600*24*365, 10000 );
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem

% T1, T2 è tspan in colonna
% Y sono le componenti posizioni e velocità
% Plot the results over 5 periods
figure(2)
Terra3d
plot3( Y1(:,1), Y1(:,2), Y1(:,3), '--k')

view([-30 20])
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit - 1 year');
axis equal;
hold on
plot3( Y2(:,1), Y2(:,2), Y2(:,3), '--r' )
legend('Unperturbed orbit','Perturbed orbit')
grid on

%% Energy 
clc
close all

tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem

r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];

epsilon1 = (vecnorm(v1,2,2)).^2/2-mu_E./vecnorm(r1,2,2); %2BP
epsilon2 = (vecnorm(v2,2,2)).^2/2-mu_E./vecnorm(r2,2,2); %Perturbed

plot(tspan,epsilon1,'-b',tspan,epsilon2,'--r','LineWidth',1.5)
grid on
xlabel('t [s]')
ylabel('\epsilon [km^2/s^2]')
title('Energy over 5 periods')
legend('2BP','J2 Perturbed')

% over 1 year
tspan = linspace( 0, 24*3600*365, 10000 );

% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem

r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];

epsilon1 = (vecnorm(v1,2,2)).^2/2-mu_E./vecnorm(r1,2,2); %2BP
epsilon2 = (vecnorm(v2,2,2)).^2/2-mu_E./vecnorm(r2,2,2); %Perturbed

plot(tspan,epsilon1,'-b',tspan,epsilon2,'--r','LineWidth',1.5)
grid on
xlabel('t [s]')
ylabel('\epsilon [km^2/s^2]')
title('Energy over 1 year')
legend('2BP','J2 Perturbed')

%% Angular momentum 
clc
close all


r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);


tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem


% Comparison between angular momentum of unperturbed vs perturbed orbit
% over 5 periods

plot(tspan,h1(:,1),'--b','LineWidth',1.5);
hold on;
plot(tspan,h1(:,2),'--r','LineWidth',1.5);
plot(tspan,h1(:,3),'--g','LineWidth',1.5);
hold on
plot(tspan,h2_norm,'--k','Linewidth',2)
grid on;
title('Momento angolare 5 periods')
xlabel('t [s]')
ylabel('h [km^2/s]')


hold on
plot(tspan,h2(:,1),'-.b','LineWidth',1.5);
hold on;
plot(tspan,h2(:,2),'-.r','LineWidth',1.5);
plot(tspan,h2(:,3),'-.g','LineWidth',1.5);
hold on
plot(tspan,h2_norm,'-.k','Linewidth',1.5)

legend('h_x, 2BP','h_y,2BP','h_z,2BP','h_x, J2 perturbed','h_y,J2 Perturbed','h_z,J2 Perturbed')

% Over 1 year period

tspan = linspace( 0, 24*3600*365, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem

r1 = [Y1(:,1) Y1(:,2) Y1(:,3)];
v1 = [Y1(:,4) Y1(:,5) Y1(:,6)];
h1 = cross(r1,v1);
h1_norm = vecnorm(h1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);
h2_norm = vecnorm(h2,2,2);

plot(tspan,h1(:,1),'--b','LineWidth',1.5);
hold on;
plot(tspan,h1(:,2),'--r','LineWidth',1.5);
plot(tspan,h1(:,3),'--g','LineWidth',1.5);
hold on
plot(tspan,h2_norm,'--k','Linewidth',2)
grid on;
title('Angular momentum over 1 year')
xlabel('t [s]')
ylabel('h [km^2/s]')


hold on
plot(tspan,h2(:,1),'-.b','LineWidth',1.5);
hold on;
plot(tspan,h2(:,2),'-.r','LineWidth',1.5);
plot(tspan,h2(:,3),'-.g','LineWidth',1.5);
hold on
plot(tspan,h2_norm,'-.k','Linewidth',1.5)

legend('h_x, 2BP','h_y,2BP','h_z,2BP','h_x, J2 perturbed','h_y,J2 Perturbed','h_z,J2 Perturbed')

%% Eccentricity vector
close all

tspan = linspace( 0, 5*T, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem

e1 = 1/mu_E*cross(v1,h1)-r1./vecnorm(r1,2,2); % eccentricity vector unperturbed orbit
e1_norm = vecnorm(e1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);

e2 = 1/mu_E*cross(v2,h2)-r2./vecnorm(r2,2,2); % eccentricity vector perturbed orbit
e2_norm = vecnorm(e2,2,2);


plot(tspan,e1(:,1),'--b',tspan,e1(:,2),'--r',tspan,e1(:,3),'--g',tspan,e1_norm,'--k','LineWidth',1.5)
grid on

hold on
plot(tspan,e2(:,1),':b',tspan,e2(:,2),':r',tspan,e2(:,3),':g',tspan,e2_norm,':k','LineWidth',1.5)

title('Eccentricity vector - 5 periods')
legend('e_x, 2BP','e_y, 2BP','e_z, 2BP','|e|','e_x, J2 perturbed','e_y, J2 perturbed','e_z, J2 perturbed','|e|, J2 perturbed')


% Over 1 year

tspan = linspace( 0, 24*3600*365, 1000 );
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration 
[ T1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options ); % solution of the unperturbed problem
[ T2, Y2 ] = ode113( @(t,y) ode_2bp_j2( t, y, mu_E,R_Earth,J2 ), tspan, y0, options ); % solution of the perturbed problem

e1 = 1/mu_E*cross(v1,h1)-r1./vecnorm(r1,2,2); % eccentricity vector unperturbed orbit
e1_norm = vecnorm(e1,2,2);

r2 = [Y2(:,1) Y2(:,2) Y2(:,3)];
v2 = [Y2(:,4) Y2(:,5) Y2(:,6)];
h2 = cross(r2,v2);

e2 = 1/mu_E*cross(v2,h2)-r2./vecnorm(r2,2,2); % eccentricity vector perturbed orbit
e2_norm = vecnorm(e2,2,2);


plot(tspan,e1(:,1),'--b',tspan,e1(:,2),'--r',tspan,e1(:,3),'--g',tspan,e1_norm,'--k','LineWidth',1.5)
grid on

hold on
plot(tspan,e2(:,1),':b',tspan,e2(:,2),':r',tspan,e2(:,3),':g',tspan,e2_norm,':k','LineWidth',1.5)

title('Eccentricity vector - 1 year')
legend('e_x, 2BP','e_y, 2BP','e_z, 2BP','|e|','e_x, J2 perturbed','e_y, J2 perturbed','e_z, J2 perturbed','|e|, J2 perturbed')

%%
close all

x = linspace(0,10,100);
y = sin(x);
plot(x,y)
grid on