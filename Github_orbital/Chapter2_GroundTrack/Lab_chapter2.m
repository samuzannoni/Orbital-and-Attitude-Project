%% Chapter 2 - Groundtrack

clc
clear all
close all

% First orbit
r0 = [-4578.219, -801.084, -7929.708]'; % Km
v0 = [0.800, -6.037, 1.385]'; % Km/s

a = 8350; % Km
e = 0.1976; % eccentricity
i = 60; % deg 
Om = 270; % deg 
om = 45; % deg 
f0 = 230; % deg: initial true anomaly

mu_E = astroConstants(13);
T = 2*pi *sqrt(a^3/mu_E);

n_orbit = 3.25;
tspan_orbit = linspace(0,n_orbit*T,1000);

%%%%%%%%%%%
y0 = [ r0; v0 ];
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_orbit, y0, options );
figure;
Terra3d
plot3(Y(:,1), Y(:,2), Y(:,3), '--');
%%%%%%%%%%%%%

% Compute the Groundtrack

w_E = deg2rad(15.04/3600); % rad/s
theta_G0 = zeros(1,length(tspan_orbit)); % Rad

t0 = zeros(1,length(tspan_orbit));

[alpha,delta,lon,lat] = groundTrack(r0,v0,theta_G0,tspan_orbit,w_E,mu_E,t0)

%% Plotting routine case 1

alpha = rad2deg(alpha);
delta = rad2deg(delta);
lat = rad2deg(lat);
lon=wrapTo180(rad2deg(lon));

% figure
% % imread("Earth_image.jpg")
% % image('Earth_image')
% % hold on
plot(lon,lat,'LineStyle','none','Marker','.');
hold on
plot(lon(1),lat(1),'O','Markersize',10)
hold on
plot(lon(1),lat(1),'oc',lon(end),lat(end),'sc','MarkerSize',10, 'LineWidth',2)
legend('Groundtrack','Start','End')
xlabel('Longitude')
ylabel('Latitude')
axis([-180 180 -90 90])

%% Case 2

clc
clear all
close all

% First orbit
r0 = [3108.128, -1040.299, -6090.022]' % Km
v0 = [5.743, 8.055, 1.555]' ; % Km/s

a = 26600; % Km
e = 0.74; % eccentricity
i = 63.4; % deg 
Om = 50; % deg 
om = 280; % deg 
f0 = 0; % deg: initial true anomaly

mu_E = astroConstants(13);
T = 2*pi *sqrt(a^3/mu_E);

n_orbit = 30;
tspan_orbit = linspace(0,n_orbit*T,1000);

%%%%%%%%%%%
y0 = [ r0; v0 ];
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_orbit, y0, options );
figure;
Terra3d
plot3(Y(:,1), Y(:,2), Y(:,3), '--');
%%%%%%%%%%%%%

% Compute the Groundtrack

w_E = deg2rad(15.04/3600); % rad/s
theta_G0 = zeros(1,length(tspan_orbit)); % Rad

t0 = zeros(1,length(tspan_orbit));

[alpha,delta,lon,lat] = groundTrack(r0,v0,theta_G0,tspan_orbit,w_E,mu_E,t0)


alpha = rad2deg(alpha);
delta = rad2deg(delta);
lat = rad2deg(lat);
lon=wrapTo180(rad2deg(lon));

figure
% imread("Earth_image.jpg")
% image('Earth_image')
% hold on
plot(lon,lat,'LineStyle','none','Marker','.');
hold on
plot(lon(1),lat(1),'O','Markersize',10)
hold on
plot(lon(1),lat(1),'oc',lon(end),lat(end),'sc','MarkerSize',10, 'LineWidth',2)
legend('Groundtrack','Start','End')
xlabel('Longitude')
ylabel('Latitude')
axis([-180 180 -90 90])


%% Case 3

clc
clear all
close all


            a = astroConstants(23)+800; % Km
            e = 0; % eccentricity
            Om = 0; % deg
            om = 40; % deg
            f0 = 0; % deg: initial true anomaly

orbit = input('Choice of orbit: ')
switch orbit
    case '1'
            r0 = [5493.312, 4609.436, 0]'; % Km
            v0 = [-4.792, 5.711, 0]' ; % Km
            i = 0; % deg
            
    case '2'
            r0 = [5493.312, 3991.889, 2304.718]'; % Km
            v0 = [-4.792, 4.946, 2.856]' ; % Km
            i = 30; % deg

    case '3'
            r0 = [5493.312, -641.510, 4564.578]' % Km
            v0 = [-4.792, -0.795, 5.656]' ; % Km
            i = 98; % deg

end


mu_E = astroConstants(13);
T = 2*pi *sqrt(a^3/mu_E);

n_orbit = 5;
tspan_orbit = linspace(0,n_orbit*T,1000);

%%%%%%%%%%%
y0 = [ r0; v0 ];
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_orbit, y0, options );
figure;
Terra3d
plot3(Y(:,1), Y(:,2), Y(:,3), '--');
%%%%%%%%%%%%%

% Compute the Groundtrack

w_E = deg2rad(15.04/3600); % rad/s
theta_G0 = zeros(1,length(tspan_orbit)); % Rad

t0 = zeros(1,length(tspan_orbit));

[alpha,delta,lon,lat] = groundTrack(r0,v0,theta_G0,tspan_orbit,w_E,mu_E,t0)


alpha = rad2deg(alpha);
delta = rad2deg(delta);
lat = rad2deg(lat);
lon=wrapTo180(rad2deg(lon));

figure
% imread("Earth_image.jpg")
% image('Earth_image')
% hold on
plot(lon,lat,'LineStyle','none','Marker','.');
hold on
plot(lon(1),lat(1),'O','Markersize',10)
hold on
plot(lon(1),lat(1),'oc',lon(end),lat(end),'sc','MarkerSize',10, 'LineWidth',2)
legend('Groundtrack','Start','End')
xlabel('Longitude')
ylabel('Latitude')
axis([-180 180 -90 90])


