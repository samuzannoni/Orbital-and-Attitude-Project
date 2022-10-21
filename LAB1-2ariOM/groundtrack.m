function [alpha,delta,lon,lat] = groundtrack( r0,v0,th_g0,n_orb)

mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]

om_E = deg2rad(15.04/3600); % 15.04 deg/h

[a, e, i, OM, om, th] = car2kep(r0,v0,mu_E);

y0 = [ r0, v0];

% Integration timespan
time_step = 10;
T = 2*pi*sqrt(a^3/mu_E);
t_span = [0:time_step:n_orb*T];

% Set integration options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[~, Y] = ode113( @(t,y) ode_2bp(t,y, mu_E), t_span, y0, options );

R = Y(:,1:3);
Rnorm = vecnorm(R')';
[alpha, delta] = car2RADec(R);

thetag_vect = wrapTo2Pi(th_g0+om_E*t_span);

lon = rad2deg(wrapTo2Pi(atan2((R(:,2)./Rnorm),(R(:,1)./Rnorm))-thetag_vect')-pi);
lat = rad2deg(delta);


%  Plot
F = figure('Name','Ground Track');
set(F, 'Units', 'Normalized', 'OuterPosition', [.15 .25 .7 .7]);
hold on;
axis equal;
set(gca,'XTick',[-180:15:180],'XTickMode','manual');
set(gca,'YTick',[-90:10:90],'YTickMode','manual');
xlim([-180,180]); ylim([-90,90]);
        
image_file = 'earth.jpg';
cdata= flip(imread(image_file));
imagesc([-180,180],[-90, 90],cdata);

plot(lon, lat, '.g');
xlabel('Longitude \lambda   [deg]')
ylabel('Latitude \phi   [deg]')

end
