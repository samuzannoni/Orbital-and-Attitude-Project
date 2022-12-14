clc
clear all
close all


addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\timeConversion\time\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\plotCelBody\');


        %%%%%%%%%%%%%%%%%%%%%% Departure time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_dep_in = date2mjd2000([2026 06 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2028 06 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_step = 24/12; % The sample is taken evry t_step hour
t_span_dep = linspace(t_dep_in,t_dep_fin,(t_dep_fin-t_dep_in)*(t_step));% time vector for departure window (it gives me the step every 2 hours)


%%%%%%%%%%%%%%%%%%%%%%% Arrival time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_arr_in = date2mjd2000([2028 06 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2034 01 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_arr = linspace(t_arr_in,t_arr_fin,(t_arr_fin-t_arr_in)*(t_step)); % time vector for departure window (it gives me the step every 2 hours)



%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Earth's ephemerides %%%%%%%%%%%%%%%%%%%%%%

DEP_eph_E = [];
for k1 = 1:length(t_span_dep)

[kep_dep_E,ksun] = uplanet(t_span_dep(k1), 3); % ephemerides of Earth at departure

a_E = kep_dep_E(1);
e = kep_dep_E(2);
i = kep_dep_E(3);
Om =  kep_dep_E(4);
om =  kep_dep_E(5);
f0 =  kep_dep_E(6);

mu = ksun; 
[r_dep_E,v_dep_E] = kep2car_mod(a_E,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_dep_E = [r_dep_E' v_dep_E']; % Matrix containing Earth's state for departure window

DEP_eph_E = [DEP_eph_E; s_dep_E]; % matrix containing all Earth's states for departure time window

end



%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Jupiter's ephemerides %%%%%%%%%%%%%%%%%%%%%%

ARR_eph_J = [];
mu = ksun;
for k2 = 1:length(t_span_arr)

[kep_arr_J,ksun] = uplanet(t_span_arr(k2), 5); % ephemerides of Jupiter at arrival window

a_J = kep_arr_J(1);
e = kep_arr_J(2);
i = kep_arr_J(3);
Om =  kep_arr_J(4);
om =  kep_arr_J(5);
f0 =  kep_arr_J(6);

[r_arr_J,v_arr_J] = kep2car_mod(a_J,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_arr_J = [r_arr_J' v_arr_J']; % Matrix containing Mars's state for arrival window


ARR_eph_J= [ARR_eph_J; s_arr_J]; % matrix containing all Mars's states for arrival time window
end

%%%%%%%%%%%%%%%%%%%%% COMPUTING OF LAMBERT ARCS AND DELTA V %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Lambert arcs between departure and arrival windows %%%%%%%%%%%%
    
MU = ksun;
DEP_eph = DEP_eph_E;
ARR_eph = ARR_eph_J;
V_inf = 9.1;
[DV_constrained,DV_unconstrained,min_V,date_dep,date_arr,t_dep_index,t_arr_index,DV1,DV2,VI_T,RI_T,RF_T] = Transfer_arc(t_span_dep,t_span_arr,MU,DEP_eph,ARR_eph,V_inf);




%%%%%% Plot the transfer arc %%%%%%

TOF = (t_span_arr(t_arr_index) - t_span_dep(t_dep_index)) *24*3600; % [s]
orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
Nrev = 0; %Zero-revolution case
Ncase = 0; %Not used for the zero-revolution case
optionsLMR = 0; %No display
RF_T = RF_T';

[a_T,P,E,ERROR,VI_T,VF,TPAR,THETA] = lambertMR(RI_T,RF_T,TOF,MU,orbitType,Nrev,Ncase,optionsLMR);


%%%%%%% Propagation of Earth's orbit %%%%%%

mu_E = astroConstants(13);
RI_E = DEP_eph_E(t_dep_index,1:3); % Earth's position vector at departure
VI_E = DEP_eph_E(t_dep_index,4:6); % Earth's velocity vector at departure

y0_E = [ RI_E'; VI_E']; % transpose because input of ode
T_E = 2*pi*sqrt(a_E^3/ksun);
tspan_E = linspace(0,T_E,1000);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ ~,Y_E ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_E, y0_E, options ); % Y is the position vector of Earth over its orbit


%%%%%%% Propagation of Jupiter's orbit %%%%%%

mu_J = astroConstants(15);
RI_J = ARR_eph_J(t_arr_index,1:3); % Jupiter's position vector at departure
VI_J = ARR_eph_J(t_arr_index,4:6); % Jupiter's position vector at departure

y0_J = [ RI_J'; VI_J']; % transpose because input of ode
T_J = 2*pi*sqrt(a_J'^3/ksun);
tspan_J = linspace(0,T_J,1000);

% Perform the integration
[ ~,Y_J ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_J, y0_J, options ); % Y is the position vector of Mars over its orbit


%%%%%%%%%%%%% Propagation of transfer arc %%%%%%%%%%%%%%%


y0_T = [ RI_T; VI_T']; % transpose because input of ode
T_T = TOF;
tspan_T = linspace(0,TOF,1000);
% Perform the integration
[ ~,Y_T ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_T, y0_T, options ); % Y is the position vector of the transfer arc


%%%%%%%%%%%% Plot of the "Pork-chop plot" %%%%%%%%%%%%%%%%

[X,Y] = meshgrid(t_span_dep,t_span_arr);


figure(1)
%surf(X,Y,DV_unconstrained);
grid on
contour(X,Y,DV_unconstrained','LineWidth',1.5)
colorbar
caxis([2 30])
%xticks([t_span_dep(1) t_span_dep(60) t_span_dep(91) t_span_dep(122) t_span_dep(152), t_span_dep(182)])
xticks([t_span_dep(1):60:t_span_dep(end)]);
xticklabels({'2026 06 01 12:00','2026 07 01 12:00','2026 08 01 12:00','2026 09 01 12:00','2026 10 01 12:00', '2026 11 01 12:00', '2026 12 01 12:00','2027 01 01 12:00','2027 02 01 12:00','2027 03 01 12:00','2027 04 01 12:00','2027 05 01 12:00','2027 06 01 12:00','2027 07 01 12:00','2027 08 01 12:00','2027 09 01 12:00','2027 10 01 12:00','2027 11 01 12:00','2027 12 01 12:00','2028 01 01 12:00','2028 02 01 12:00','2028 03 01 12:00','2028 04 01 12:00','2028 05 01 12:00','2028 06 01 12:00'})
set(gca,'XTickLabelRotation',45)

%yticks([t_span_arr(1) t_span_arr(360) t_span_arr(732) t_span_arr(1092) t_span_arr(1464) t_span_arr(1836) t_span_arr(end)])
yticks([t_span_arr(1):60:t_span_arr(end)]);
yticklabels({'2026 06 01 12:00','2026 07 01 12:00','2026 08 01 12:00','2026 09 01 12:00','2026 10 01 12:00', '2026 11 01 12:00', '2026 12 01 12:00','2027 01 01 12:00','2027 02 01 12:00','2027 03 01 12:00','2027 04 01 12:00','2027 05 01 12:00','2027 06 01 12:00','2027 07 01 12:00','2027 08 01 12:00','2027 09 01 12:00','2027 10 01 12:00','2027 11 01 12:00','2027 12 01 12:00','2028 01 01 12:00','2028 02 01 12:00','2028 03 01 12:00','2028 04 01 12:00','2028 05 01 12:00','2028 06 01 12:00','2028 07 01 12:00','2028 08 01 12:00','2028 09 01 12:00','2028 10 01 12:00','2028 11 01 12:00', '2028 12 01 12:00', '2029 01 01 12:00','2029 02 01 12:00','2029 03 01 12:00','2029 04 01 12:00','2029 05 01 12:00','2029 06 01 12:00','2029 07 01 12:00','2029 08 01 12:00','2029 09 01 12:00','2029 10 01 12:00','2029 11 01 12:00','2029 12 01 12:00','2030 01 01 12:00','2030 02 01 12:00','2030 03 01 12:00','2030 04 01 12:00','2030 05 01 12:00','2030 06 01 12:00','2030 07 01 12:00','2030 08 01 12:00','2030 09 01 12:00','2030 10 01 12:00','2030 11 01 12:00','2030 12 01 12:00', '2031 01 01 12:00', '2031 02 01 12:00','2031 03 01 12:00','2031 04 01 12:00','2031 05 01 12:00','2031 06 01 12:00','2031 07 01 12:00','2031 08 01 12:00','2031 09 01 12:00','2031 10 01 12:00','20321 11 01 12:00','2031 12 01 12:00','2027 11 01 12:00','2027 12 01 12:00','2028 01 01 12:00','2028 02 01 12:00','2028 03 01 12:00','2028 04 01 12:00','2028 05 01 12:00','2028 06 01 12:00'})
set(gca,'YTickLabelRotation',45)

%%%%%%%%% Plot the transfer arc %%%%%%

figure(2) 
plotCelBody(0,50e6)
hold on
plot3( Y_E(:,1),Y_E(:,2), Y_E(:,3),'LineWidth',1.5); % Plot of initial orbit (Earth)
plot3(Y_J(:,1),Y_J(:,2),Y_J(:,3),'LineWidth',1.5); % Plot of final orbit (Mars)

plotCelBody(3,10e6,[RI_E(1),RI_E(2),RI_E(3)]) % plot of the Earth (not scaled)
plotCelBody(5,35e6,[RI_J(1),RI_J(2),RI_J(3)]) % plot of Jupiter (not scaled)
% plot3(RI_E(1),RI_E(2),RI_E(3),'ob','MarkerFaceColor','b','MarkerSize',10); % Plot of P1
% plot3(RI_M(1),RI_M(2),RI_M(3),'or','MarkerFaceColor','r','MarkerSize',10) % Plot of P2
axis equal
grid on
comet3(Y_T(:,1),Y_T(:,2),Y_T(:,3)) % Animation of transfer orbit
legend('Sun','Earth''s orbit','Jupiter''s orbit','Earth''s position at departure','Jupiter''s position at arrival','Transfer arc')
title('Transfer trajectory from Earth to Mars')
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')