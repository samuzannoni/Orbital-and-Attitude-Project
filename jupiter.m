%% JUPITER(5)
% Departure time window 
t_dep_in = date2mjd2000([2026 06 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2028 06 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00
n1 =(t_dep_fin-t_dep_in+1); % step time span

t_span_dep = linspace(t_dep_in,t_dep_fin,n1); % time vector for departure window (be careful it's in days!!!)

% Arrival time window 
t_arr_in = date2mjd2000([2028 06 01 12 00 00]); %initial arrival time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2034 01 01 12 00 00]); %final arrival time given in number of days since 01 01 2000 at 12:00
n2= (t_arr_fin-t_arr_in+1);

t_span_arr = linspace(t_arr_in,t_arr_fin,n2);% time vector for arrival window 

mu_S = astroConstants(4); %Gravitational Constant of the Sun

V_inf= 9.1;

[DV_tot,DV1,min_DV_tot,dv1_verificata,data1,data2,ROW,COL] = Deltavtotconstrain(3,5,t_span_dep,t_span_arr,mu_S,n1,n2,V_inf);
%
%[DV_tot2,min_DV_tot2,date12,date22] = Deltavtot(3,5,t_span_dep,t_span_arr,mu_S,n1,n2);
%
figure
% subplot(1,2,1)
[X1, X2] = meshgrid(t_span_dep, t_span_arr);
contour(X1,X2,DV_tot',[14.6 14.8],'LineWidth', 1) % 'ShowText', 'on');
hold on;
plot(t_span_dep(ROW),t_span_arr(COL),'ob','MarkerFaceColor','b','MarkerSize',6)
xlabel('Departure')
ylabel('Arrival')
X3 = colorbar;
caxis([14 20]);
grid on;