clc;
close all;
clear all;

addpath ('C:\Users\arilu\Documents\LAB-OM\LAB1-2\time')
n = 1300;
% Departure time window 
t_dep_in = date2mjd2000([2003 04 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2003 08 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_dep = linspace(t_dep_in,t_dep_fin,n); % time vector for departure window (be careful it's in days!!!)

% Arrival time window 
t_arr_in = date2mjd2000([2003 09 01 12 00 00]); %initial arrival time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2004 03 01 12 00 00]); %final arrival time given in number of days since 01 01 2000 at 12:00

t_span_arr = linspace(t_arr_in,t_arr_fin,n);% time vector for arrival window (be careful it's in days!!!)

MU = astroConstants(4); 

% Compute the total cost of the manoeuvre
[DV_tot,min_DV_tot,date1,date2] = vtot(3,4,t_span_dep,t_span_arr,MU,n);

%% PorkChop Plot
figure

% subplot(1,2,1)
[X1, X2] = meshgrid(t_span_arr, t_span_dep);
% surf(X1, X2, DV_tot, 'LineStyle', 'None');
% xlabel('TOF')
% ylabel('tspandep')
% zlabel('DVtot')
% zlim([min_DV_tot 15])
% X3 = colorbar;
% 
% subplot(1,2,2)
contour(X1,X2,DV_tot,100, 'LineWidth', 1) % 'ShowText', 'on');
xlabel('TOF')
ylabel('tspandep')
X3 = colorbar;
caxis([5 10]);




