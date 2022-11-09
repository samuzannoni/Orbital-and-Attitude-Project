t_epoch = horizonsresults(:,1);
      e = horizonsresults(:,2);
      i = horizonsresults(:,3);
     OM = horizonsresults(:,4);
     om = horizonsresults(:,5);
     TA = horizonsresults(:,6); % true anomaly
      A = horizonsresults(:,7); % semi-major axis

clc
close all

addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\Image_celestial_bodies\');


%% punto 4
mu_E = astroConstants(13);
w_E = deg2rad(15.04/3600);

rr = [];
vv = [];


for j = 1:length(t_epoch)

    [r, v] = kep2car(e(j),OM(j),i(j),om(j),TA(j),A(j));
    
    rr = [rr; r'];
    vv = [vv; v'];

end


theta_G0 = zeros(1,length(t_epoch));
tspan_orbit = linspace(0, 24*3600, length(t_epoch));
t0 = zeros(1,length(t_epoch));

[lon,lat] = groundTrackMOD(rr,vv,theta_G0,tspan_orbit,w_E,mu_E,t0);


lat = rad2deg(lat);
lon = wrapTo180(rad2deg(lon));

lon_plot = lon;
lat_plot = lat;

for z=1:length(lon)-1
    if (lon_plot(z)>0 && lon_plot(z+1)<0)
        lon_plot(z) = NaN;
        lat_plot(z) = NaN;
    end
end


close all
figure(1);
E_image = imread("Earth.jpg");
plot(lon_plot, lat_plot,'-r');

hold on
plot(lon_plot(1),lat_plot(1),'oc',lon_plot(end),lat_plot(end),'sc','MarkerSize',10, 'LineWidth',2, 'MarkerEdgeColor','r');
axis([-180 180 -90 90]);
Im = image(xlim, -ylim, E_image);

xlabel('Longitude');
ylabel('Latitude');
axis([-180 180 -90 90]);

uistack(Im,'bottom');


%% punto 5 

t = horizonsresults(:,1);
e = horizonsresults(:,2);
i = horizonsresults(:,3);
OM = horizonsresults(:,4);
om = horizonsresults(:,5);
TA = horizonsresults(:,6); % true anomaly
A = horizonsresults(:,7); % semi-major axis

mu_E = astroConstants(13);
w_E = deg2rad(15.04/3600);


tspan_orbit = [];
t = t*24*3600; %s, converte in secondi i julian days

for k = 1:length(t)

    tspan_orbit = [tspan_orbit; t(k) - t(1)];

end

tspan_orbit = tspan_orbit';

theta_G0=zeros(1, length(tspan_orbit));
t0=zeros(1, length(tspan_orbit));

[r0, v0] = kep2car(e(1), OM(1), i(1), om(1), TA(1), A(1));
[alpha,delta,lon,lat] = groundTrack(r0,v0,theta_G0,tspan_orbit,w_E,mu_E,t0);

lat=rad2deg(lat);
lon=wrapTo180(rad2deg(lon));


lon_plot = lon;
lat_plot = lat;

for z = 1:length(lon)-1
    if (lon_plot(z)>0 && lon_plot(z+1)<0)
        lon_plot(z) = NaN;
        lat_plot(z) = NaN;
    end
end

plot(lon_plot, lat_plot,'-g');
plot(lon_plot(end),lat_plot(end),'sc','MarkerSize',10, 'LineWidth',2,'MarkerEdgeColor','g');

legend('Groundtrack from ephemerides', 'Start_1', 'End_1','Groundtrack of propagated orbit','End_2');