%%5.
clc
close all

addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\');

mu_E = astroConstants(13);
w_E = deg2rad(15.04/3600);

e = horizonsresults(:,2);
i = horizonsresults(:,3);
OM = horizonsresults(:,4);
om = horizonsresults(:,5);
TA = horizonsresults(:,6); % true anomaly
A = horizonsresults(:,7); % semi-major axis

mu_E = astroConstants(13);
w_E = deg2rad(15.04/3600);

n_orbit=15.5;
T=2*pi*sqrt(A(1)^3/mu_E);

tspan_orbit=linspace(0, n_orbit*T, 24*60);
theta_G0=zeros(1, length(tspan_orbit));
t0=zeros(1, length(tspan_orbit));

[r0, v0] = kep2car(e(1), OM(1), i(1), OM(1), TA(1), A(1));
[alpha,delta,lon,lat] = groundTrack(r0,v0,theta_G0,tspan_orbit,w_E,mu_E,t0);

lat=rad2deg(lat);
lon=wrapTo180(rad2deg(lon));

figure;
plot(lon, lat, 'LineStyle', 'none', 'Marker', '.');
hold on
plot(lon(1),lat(1),'oc',lon(end),lat(end),'sc','MarkerSize',10, 'LineWidth',2);
legend('Groundtrack', 'Start', 'End');
xlabel('Longitude');
ylabel('Latitude');
axis([-180 180 -90 90]);