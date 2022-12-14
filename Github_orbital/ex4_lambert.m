%%%%%% Excercise 4 of Chapter 3 %%%%%%%


clc
clear all
close all 

addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\timeConversion\time\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\plotCelBody\');

%%%% Same excercise of ex3_Lambert but with the launcher constraint %%%%%%%
%%%% Problem with pork-chop plot because DV_constrained is not a matrix but
%%%% a vector
%%%% The dimensions are not scaled (has to be fixed)


%%% DATA %%%%
R_Sun = 695000*50; % [Km]
R_Earth = 6378; % [Km]
R_Mars = 3389; % [Km]
Transfer = input('Where do you want to go? :) : ')

switch Transfer


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            %Case 1: Mercury
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Mercury'





%%%%%%%%%%%%%%%%%%%%%% Departure time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_dep_in = date2mjd2000([2023 01 11 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2025 01 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_dep = linspace(t_dep_in,t_dep_fin,1000);% time vector for departure window (it gives me the step every 2 hours)


%%%%%%%%%%%%%%%%%%%%%%% Arrival time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_arr_in = date2mjd2000([2024 04 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2025 03 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_arr = linspace(t_arr_in,t_arr_fin,1000); % time vector for departure window (it gives me the step every 2 hours)
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



%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Mercury's ephemerides %%%%%%%%%%%%%%%%%%%%%%

ARR_eph_Me = [];
mu = ksun;
for k2 = 1:length(t_span_arr)

[kep_arr_Me,ksun] = uplanet(t_span_arr(k2), 1); % ephemerides of Mars at arrival window

a_Me = kep_arr_Me(1);
e = kep_arr_Me(2);
i = kep_arr_Me(3);
Om =  kep_arr_Me(4);
om =  kep_arr_Me(5);
f0 =  kep_arr_Me(6);

[r_arr_Me,v_arr_Me] = kep2car_mod(a_Me,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_arr_Me = [r_arr_Me' v_arr_Me']; % Matrix containing Mercury's state for arrival window


ARR_eph_Me = [ARR_eph_Me; s_arr_Me]; % matrix containing all Mercury's states for arrival time window
end

%%%%%%%%%%%%%%%%%%%%% COMPUTING OF LAMBERT ARCS AND DELTA V %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Lambert arcs between departure and arrival windows %%%%%%%%%%%%
    
MU = ksun;
DEP_eph = DEP_eph_E;
ARR_eph = ARR_eph_Me;
V_inf = 7; % [Km/s], Launcher constraint

[DV_constrained,DV_unconstrained,min_V,date_dep,date_arr,t_dep_index,t_arr_index,DV1,DV2,VI_T,RI_T,RF_T] = Transfer_arc(t_span_dep,t_span_arr,MU,DEP_eph,ARR_eph,V_inf);

%%%%%%%%%%%% Plot of the "Pork-chop plot" %%%%%%%%%%%%%%%%

[X,Y] = meshgrid(t_span_dep,t_span_arr);
close all

% figure(1)
% surf(X,Y,DV_constrained);
% colorbar;
% grid on
% hold on
% contour(X,Y,DV_constrained,100,'LineWidth',1.5)
%xticks([t_span_dep(1) t_span_dep(245) t_span_dep(499) t_span_dep(745) t_span_dep(1000)])
%xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'})
%set(gca,'XTickLabelRotation',45)

% yticks([t_span_arr(1) t_span_arr(165) t_span_arr(335) t_span_arr(500) t_span_arr(670) t_span_arr(841) t_span_arr(1000)])
% yticklabels({'2003 Sep 01','2003 Oct 01','2003 Nov 01','2003 Dec 01','2004 Jan 01','2004 Feb 01','2004 Mar 01'})
% set(gca,'YTickLabelRotation',45)
%colorbar


%%%%%% Plot the transfer arc %%%%%%



TOF = (t_span_arr(t_arr_index) - t_span_dep(t_dep_index)) *24*3600; % [s]
orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
Nrev = 0; %Zero-revolution case
Ncase = 0; %Not used for the zero-revolution case
optionsLMR = 0; %No display

[a_T,P,E,ERROR,VI_T,VF,TPAR,THETA] = lambertMR(RI_T,RF_T,TOF,MU,orbitType,Nrev,Ncase,optionsLMR);


%%%%%%% Propagation of Earth's orbit %%%%%%

mu_E = astroConstants(13);
RI_E = DEP_eph_E(t_dep_index,1:3);
VI_E = DEP_eph_E(t_dep_index,4:6);

y0_E = [ RI_E'; VI_E']; % transpose because input of ode
T_E = 2*pi*sqrt(a_E^3/ksun);
tspan_E = linspace(0,T_E,1000);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ ~,Y_E ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_E, y0_E, options ); % Y is the position vector of Earth over its orbit


%%%%%%% Propagation of Mercurys's orbit %%%%%%

mu_Me = astroConstants(11); % Mercury's gravity constant
RI_Me = ARR_eph_Me(t_arr_index,1:3);
VI_Me = ARR_eph_Me(t_arr_index,4:6);

y0_Me = [ RI_Me'; VI_Me']; % transpose because input of ode
T_Me = 2*pi*sqrt(a_Me^3/ksun);
tspan_Me = linspace(0,T_Me,1000);

% Perform the integration
[ ~,Y_Me ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_Me, y0_Me, options ); % Y is the position vector of Mars over its orbit

%%%%%%%%%%%%% Propagation of transfer arc %%%%%%%%%%%%%%%




y0_T = [ RI_T; VI_T']; % transpose because input of ode
T_T = TOF;
tspan_T = linspace(0,T_T,1000);
% Perform the integration
[ ~,Y_T ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_T, y0_T, options ); % Y is the position vector of the transfer arc

%% %%%%%% Plot the transfer arc %%%%%%

figure(3) 
plotCelBody(0,18e6)
hold on
plot3( Y_E(:,1),Y_E(:,2), Y_E(:,3),'LineWidth',1.5); % Plot of initial orbit (Earth)
hold on
plot3(Y_Me(:,1),Y_Me(:,2),Y_Me(:,3),'LineWidth',1.5); % Plot of final orbit (Mars)

plotCelBody(3,6378e3,[RI_E(1),RI_E(2),RI_E(3)])
plotCelBody(1,2439e3,[RI_Me(1),RI_Me(2),RI_Me(3)])
% plot3(RI_E(1),RI_E(2),RI_E(3),'ob','MarkerFaceColor','b','MarkerSize',10); % Plot of P1
% plot3(RI_Me(1),RI_Me(2),RI_Me(3),'oc','MarkerFaceColor','c','MarkerSize',10) % Plot of P2
axis equal
grid on
comet3(Y_T(:,1),Y_T(:,2),Y_T(:,3)) % Animation of transfer orbit
legend('Sun','Earth''s orbit','Mercurys''s orbit','Earth''s position at departure','Mercurys''s position at arrival','Transfer arc')
title('Transfer trajectory from Earth to Mercury')
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            %Case 2: Venus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    case 'Venus'

        %%%%%%%%%%%%%%%%%%%%%% Departure time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_dep_in = date2mjd2000([2024 06 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2026 11 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_dep = linspace(t_dep_in,t_dep_fin,1000);% time vector for departure window (it gives me the step every 2 hours)


%%%%%%%%%%%%%%%%%%%%%%% Arrival time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_arr_in = date2mjd2000([2024 12 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2027 06 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_arr = linspace(t_arr_in,t_arr_fin,1000); % time vector for departure window (it gives me the step every 2 hours)


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



%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Venus's ephemerides %%%%%%%%%%%%%%%%%%%%%%

ARR_eph_Ve = [];
mu = ksun;
for k2 = 1:length(t_span_arr)

[kep_arr_Ve,ksun] = uplanet(t_span_arr(k2), 2); % ephemerides of Mars at arrival window

a_Ve = kep_arr_Ve(1);
e = kep_arr_Ve(2);
i = kep_arr_Ve(3);
Om =  kep_arr_Ve(4);
om =  kep_arr_Ve(5);
f0 =  kep_arr_Ve(6);

[r_arr_Ve,v_arr_Ve] = kep2car_mod(a_Ve,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_arr_Ve = [r_arr_Ve' v_arr_Ve']; % Matrix containing Mercury's state for arrival window


ARR_eph_Ve = [ARR_eph_Ve; s_arr_Ve]; % matrix containing all Mercury's states for arrival time window
end

%%%%%%%%%%%%%%%%%%%%% COMPUTING OF LAMBERT ARCS AND DELTA V %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Lambert arcs between departure and arrival windows %%%%%%%%%%%%
    
MU = ksun;
DEP_eph = DEP_eph_E;
ARR_eph = ARR_eph_Ve;
V_inf = 3; % [Km/s], Launcher constraint
[DV_constrained,DV_unconstrained,min_V,date_dep,date_arr,t_dep_index,t_arr_index,DV1,DV2,VI_T,RI_T,RF_T] = Transfer_arc(t_span_dep,t_span_arr,MU,DEP_eph,ARR_eph,V_inf);

%%%%%%%%%%%% Plot of the "Pork-chop plot" %%%%%%%%%%%%%%%%

[X,Y] = meshgrid(t_span_dep,t_span_arr);
close all

% figure(1)
% surf(X,Y,DV_constrained);
% colorbar;
% grid on
% hold on
% contour(X,Y,DV_constrained,100,'LineWidth',1.5);
% xticks([t_span_dep(1) t_span_dep(245) t_span_dep(499) t_span_dep(745) t_span_dep(1000)])
% xticklabels({'2003 Apr 01','2003 May 01','2003 Jun 01','2003 Jul 01','2003 Aug 01'})
% set(gca,'XTickLabelRotation',45)

% yticks([t_span_arr(1) t_span_arr(165) t_span_arr(335) t_span_arr(500) t_span_arr(670) t_span_arr(841) t_span_arr(1000)])
% yticklabels({'2003 Sep 01','2003 Oct 01','2003 Nov 01','2003 Dec 01','2004 Jan 01','2004 Feb 01','2004 Mar 01'})
% set(gca,'YTickLabelRotation',45)
%colorbar


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
RI_E = DEP_eph_E(t_dep_index,1:3);
VI_E = DEP_eph_E(t_dep_index,4:6);

y0_E = [ RI_E'; VI_E']; % transpose because input of ode
T_E = 2*pi*sqrt(a_E^3/ksun);
tspan_E = linspace(0,T_E,1000);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ ~,Y_E ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_E, y0_E, options ); % Y is the position vector of Earth over its orbit


%%%%%%% Propagation of Venus's orbit %%%%%%

mu_Ve = astroConstants(12); % Mercury's gravity constant
RI_Ve = ARR_eph_Ve(t_arr_index,1:3);
VI_Ve = ARR_eph_Ve(t_arr_index,4:6);

y0_Ve = [ RI_Ve'; VI_Ve']; % transpose because input of ode
T_Ve = 2*pi*sqrt(a_Ve^3/ksun);
tspan_Ve = linspace(0,T_Ve,1000);

% Perform the integration
[ ~,Y_Ve ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_Ve, y0_Ve, options ); % Y is the position vector of Mars over its orbit

%%%%%%%%%%%%% Propagation of transfer arc %%%%%%%%%%%%%%%



y0_T = [ RI_T; VI_T']; % transpose because input of ode
T_T = TOF;
tspan_T = linspace(0,TOF,1000);
% Perform the integration
[ ~,Y_T ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_T, y0_T, options ); % Y is the position vector of the transfer arc

%% %%%%%% Plot the transfer arc %%%%%%

figure(3) 
plotCelBody(0,18e6);
hold on
plot3( Y_E(:,1),Y_E(:,2), Y_E(:,3),'LineWidth',1.5); % Plot of initial orbit (Earth)
hold on
plot3(Y_Ve(:,1),Y_Ve(:,2),Y_Ve(:,3),'LineWidth',1.5); % Plot of final orbit (Mars)

plotCelBody(3,6378e3,[RI_E(1),RI_E(2),RI_E(3)]); % Earth's plot
plotCelBody(2,2439e3,[RI_Ve(1),RI_Ve(2),RI_Ve(3)]); % Venus's plot
% plot3(RI_E(1),RI_E(2),RI_E(3),'ob','MarkerFaceColor','b','MarkerSize',10); % Plot of P1
% plot3(RI_Me(1),RI_Me(2),RI_Me(3),'oc','MarkerFaceColor','c','MarkerSize',10) % Plot of P2
axis equal
grid on
comet3(Y_T(:,1),Y_T(:,2),Y_T(:,3)) % Animation of transfer orbit
legend('Sun','Earth''s orbit','Venus''s orbit','Earth''s position at departure','Venus''s position at arrival','Transfer arc')
title('Transfer trajectory from Earth to Venus')
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')
 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            %Case 3: Mars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    case 'Mars'

        %%%%%%%%%%%%%%%%%%%%%% Departure time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_dep_in = date2mjd2000([2025 08 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2031 01 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_dep = linspace(t_dep_in,t_dep_fin,(t_dep_fin-t_dep_in)*12);% time vector for departure window (it gives me the step every 2 hours)


%%%%%%%%%%%%%%%%%%%%%%% Arrival time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_arr_in = date2mjd2000([2026 11 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2032 01 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_arr = linspace(t_arr_in,t_arr_fin,(t_arr_fin-t_arr_in)*12); % time vector for departure window (it gives me the step every 2 hours)



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



%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Mars's ephemerides %%%%%%%%%%%%%%%%%%%%%%

ARR_eph_M = [];
mu = ksun;
for k2 = 1:length(t_span_arr)

[kep_arr_M,ksun] = uplanet(t_span_arr(k2), 4); % ephemerides of Mars at arrival window

a_M = kep_arr_M(1);
e = kep_arr_M(2);
i = kep_arr_M(3);
Om =  kep_arr_M(4);
om =  kep_arr_M(5);
f0 =  kep_arr_M(6);

[r_arr_M,v_arr_M] = kep2car_mod(a_M,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_arr_M = [r_arr_M' v_arr_M']; % Matrix containing Mars's state for arrival window


ARR_eph_M = [ARR_eph_M; s_arr_M]; % matrix containing all Mars's states for arrival time window
end

%%
%%%%%%%%%%%%%%%%%%%%% COMPUTING OF LAMBERT ARCS AND DELTA V %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Lambert arcs between departure and arrival windows %%%%%%%%%%%%
    
MU = ksun;
DEP_eph = DEP_eph_E;
ARR_eph = ARR_eph_M;
V_inf = 3.5;
[DV_constrained,DV_unconstrained,min_V,date_dep,date_arr,t_dep_index,t_arr_index,DV1,DV2,VI_T,RI_T,RF_T] = Transfer_arc(t_span_dep,t_span_arr,MU,DEP_eph,ARR_eph,V_inf);






%%%%%% Compute the transfer arc %%%%%%



TOF = (t_span_arr(t_arr_index) - t_span_dep(t_dep_index)) *24*3600; % [s]
orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
Nrev = 0; %Zero-revolution case
Ncase = 0; %Not used for the zero-revolution case
optionsLMR = 0; %No display
RF_T = RF_T'; % transpose because input of Lambert

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





%%%%%%% Propagation of Mars's orbit %%%%%%

mu_M = astroConstants(14);
RI_M = ARR_eph_M(t_arr_index,1:3); % Mars's position vector at departure
VI_M = ARR_eph_M(t_arr_index,4:6); % Mars's position vector at departure

y0_M = [ RI_M'; VI_M']; % transpose because input of ode
T_M = 2*pi*sqrt(a_M^3/ksun);
tspan_M = linspace(0,T_M,1000);

% Perform the integration
[ ~,Y_M ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_M, y0_M, options ); % Y is the position vector of Mars over its orbit



%%%%%%%%%%%%% Propagation of transfer arc %%%%%%%%%%%%%%%


y0_T = [ RI_T; VI_T']; % transpose because input of ode
T_T = TOF;
tspan_T = linspace(0,T_T,1000);
% Perform the integration
[ ~,Y_T ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_T, y0_T, options ); % Y is the position vector of the transfer arc




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT OF THE PORK CHOP PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [X,Y] = meshgrid(t_span_dep,t_span_arr);
% close all
% 
% figure(1)
% surf(X,Y,DV_unconstrained);
% colorbar;
% 
% figure(1)
% grid on
% contour(X,Y,DV',100,'LineWidth',1.5)
% colorbar
% caxis([2 10])
% xticks([t_span_dep(1) t_span_dep(360) t_span_dep(744) t_span_dep(1104) t_span_dep(end)])
% xticklabels({'2003 04 01 12:00','2003 05 01 12:00','2003 06 01 12:00','2003 07 01 12:00','2003 08 01 12:00'})
% set(gca,'XTickLabelRotation',45)
% 
% yticks([t_span_arr(1) t_span_arr(360) t_span_arr(732) t_span_arr(1092) t_span_arr(1464) t_span_arr(1836) t_span_arr(end)])
% yticklabels({'2003 09 01 12:00','2003 10 01 12:00','2003 11 01 00','2003 Dec 01','2004 Jan 01','2004 Feb 01','2004 Mar 01'})
% set(gca,'YTickLabelRotation',45)
%%%%%%%%% Plot the transfer arc %%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2) 
plotCelBody(0,R_Sun)
hold on
plot3( Y_E(:,1),Y_E(:,2), Y_E(:,3),'LineWidth',1.5); % Plot of initial orbit (Earth)
plot3(Y_M(:,1),Y_M(:,2),Y_M(:,3),'LineWidth',1.5); % Plot of final orbit (Mars)

plotCelBody(3,R_Earth*1000,[RI_E(1),RI_E(2),RI_E(3)]) % plot of the Earth (not scaled)
plotCelBody(4,R_Mars*2000',[RI_M(1),RI_M(2),RI_M(3)]) % plot of Mars (not scaled)
% plot3(RI_E(1),RI_E(2),RI_E(3),'ob','MarkerFaceColor','b','MarkerSize',10); % Plot of P1
% plot3(RI_M(1),RI_M(2),RI_M(3),'or','MarkerFaceColor','r','MarkerSize',10) % Plot of P2
axis equal
grid on
comet3(Y_T(:,1),Y_T(:,2),Y_T(:,3)) % Animation of transfer orbit
legend('Sun','Earth''s orbit','Mars''s orbit','Earth''s position at departure','Mars''s position at arrival','Transfer arc')
title('Transfer trajectory from Earth to Mars')
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            %Case 4: Jupiter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



case 'Jupiter'

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
close all

figure(1)
%surf(X,Y,DV_unconstrained);
grid on
contour(X,Y,DV_unconstrained','LineWidth',1.5)
colorbar
caxis([1 60']);
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








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            %Case 5: Saturn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    case 'Saturn'

        %%%%%%%%%%%%%%%%%%%%%% Departure time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_dep_in = date2mjd2000([2027 09 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2029 10 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_dep = linspace(t_dep_in,t_dep_fin,1000);% time vector for departure window (it gives me the step every 2 hours)


%%%%%%%%%%%%%%%%%%%%%%% Arrival time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_arr_in = date2mjd2000([2030 04 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2036 03 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_arr = linspace(t_arr_in,t_arr_fin,1000); % time vector for departure window (it gives me the step every 2 hours)



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



%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Saturn's ephemerides %%%%%%%%%%%%%%%%%%%%%%

ARR_eph_S = [];
mu = ksun;
for k2 = 1:length(t_span_arr)

[kep_arr_S,ksun] = uplanet(t_span_arr(k2), 6); % ephemerides of Saturn at arrival window

a_S = kep_arr_S(1);
e = kep_arr_S(2);
i = kep_arr_S(3);
Om =  kep_arr_S(4);
om =  kep_arr_S(5);
f0 =  kep_arr_S(6);

[r_arr_S,v_arr_S] = kep2car_mod(a_S,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_arr_S = [r_arr_S' v_arr_S']; % Matrix containing Saturn's state for arrival window


ARR_eph_S= [ARR_eph_S; s_arr_S]; % matrix containing all Saturn's states for arrival time window
end

%%%%%%%%%%%%%%%%%%%%% COMPUTING OF LAMBERT ARCS AND DELTA V %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Lambert arcs between departure and arrival windows %%%%%%%%%%%%
    
MU = ksun;
DEP_eph = DEP_eph_E;
ARR_eph = ARR_eph_S;
V_inf = 11.5;
[DV_constrained,DV_unconstrained,min_V,date_dep,date_arr,t_dep_index,t_arr_index,DV1,DV2,VI_T,RI_T,RF_T] = Transfer_arc(t_span_dep,t_span_arr,MU,DEP_eph,ARR_eph,V_inf);


%%%%%%%%%%%% Plot of the "Pork-chop plot" %%%%%%%%%%%%%%%%

[X,Y] = meshgrid(t_span_dep,t_span_arr);
close all

% figure(1)
% surf(X,Y,DV);
% colorbar;

% figure(1)
% grid on
% contour(X,Y,DV',100,'LineWidth',1.5)
% colorbar
% caxis([5 10])
% xticks([t_span_dep(1) t_span_dep(360) t_span_dep(744) t_span_dep(1104) t_span_dep(end)])
% xticklabels({'2003 04 01 12:00','2003 05 01 12:00','2003 06 01 12:00','2003 07 01 12:00','2003 08 01 12:00'})
% set(gca,'XTickLabelRotation',45)
% 
% yticks([t_span_arr(1) t_span_arr(360) t_span_arr(732) t_span_arr(1092) t_span_arr(1464) t_span_arr(1836) t_span_arr(end)])
% yticklabels({'2003 09 01 12:00','2003 10 01 12:00','2003 11 01 00','2003 Dec 01','2004 Jan 01','2004 Feb 01','2004 Mar 01'})
% set(gca,'YTickLabelRotation',45)


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


%%%%%%% Propagation of Saturn's orbit %%%%%%

mu_S = astroConstants(16);
RI_S = ARR_eph_S(t_arr_index,1:3); % Jupiter's position vector at departure
VI_S = ARR_eph_S(t_arr_index,4:6); % Jupiter's position vector at departure

y0_S = [ RI_S'; VI_S']; % transpose because input of ode
T_S = 2*pi*sqrt(a_S'^3/ksun);
tspan_S = linspace(0,T_S,1000);

% Perform the integration
[ ~,Y_S ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_S, y0_S, options ); % Y is the position vector of Mars over its orbit

%%%%%%%%%%%%% Propagation of transfer arc %%%%%%%%%%%%%%%


y0_T = [ RI_T; VI_T']; % transpose because input of ode
T_T = TOF;
tspan_T = linspace(0,T_T,1000);
% Perform the integration
[ ~,Y_T ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_T, y0_T, options ); % Y is the position vector of the transfer arc

%% %%%%%% Plot the transfer arc %%%%%%

figure(2) 
plotCelBody(0,20e6)
hold on
plot3( Y_E(:,1),Y_E(:,2), Y_E(:,3),'LineWidth',1.5); % Plot of initial orbit (Earth)
plot3(Y_S(:,1),Y_S(:,2),Y_S(:,3),'LineWidth',1.5); % Plot of final orbit (Mars)

plotCelBody(3,8e6,[RI_E(1),RI_E(2),RI_E(3)]); % plot of the Earth (not scaled)
plotCelBody(6,18e6,[RI_S(1),RI_S(2),RI_S(3)]); % plot of Jupiter (not scaled)
% plot3(RI_E(1),RI_E(2),RI_E(3),'ob','MarkerFaceColor','b','MarkerSize',10); % Plot of P1
% plot3(RI_M(1),RI_M(2),RI_M(3),'or','MarkerFaceColor','r','MarkerSize',10) % Plot of P2
axis equal
grid on
comet3(Y_T(:,1),Y_T(:,2),Y_T(:,3)) % Animation of transfer orbit
legend('Sun','Earth''s orbit','Saturn''s orbit','Earth''s position at departure','Saturn''s position at arrival','Transfer arc')
title('Transfer trajectory from Earth to Mars')
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            %Case 6: Uranus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


case 'Uranus'

t_dep_in = date2mjd2000([2027 01 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2029 01 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_dep = linspace(t_dep_in,t_dep_fin,1000);% time vector for departure window (it gives me the step every 2 hours)


%%%%%%%%%%%%%%%%%%%%%%% Arrival time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_arr_in = date2mjd2000([2031 04 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2045 03 12 01 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_arr = linspace(t_arr_in,t_arr_fin,1000); % time vector for departure window (it gives me the step every 2 hours)






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


%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Uranus's ephemerides %%%%%%%%%%%%%%%%%%%%%%

ARR_eph_U = [];
mu = ksun;
for k2 = 1:length(t_span_arr)

[kep_arr_U,ksun] = uplanet(t_span_arr(k2), 7); % ephemerides of Saturn at arrival window

a_U = kep_arr_U(1);
e = kep_arr_U(2);
i = kep_arr_U(3);
Om =  kep_arr_U(4);
om =  kep_arr_U(5);
f0 =  kep_arr_U(6);

[r_arr_U,v_arr_U] = kep2car_mod(a_U,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_arr_U = [r_arr_U' v_arr_U']; % Matrix containing Saturn's state for arrival window


ARR_eph_U = [ARR_eph_U; s_arr_U]; % matrix containing all Saturn's states for arrival time window


end

MU = ksun;
DEP_eph = DEP_eph_E;
ARR_eph = ARR_eph_U;
V_inf = 12.1;

[DV_constrained,DV_unconstrained,min_V,date_dep,date_arr,t_dep_index,t_arr_index,DV1,DV2,VI_T,RI_T,RF_T] = Transfer_arc(t_span_dep,t_span_arr,MU,DEP_eph,ARR_eph,V_inf);

%%%%%%%%%%%% Plot of the "Pork-chop plot" %%%%%%%%%%%%%%%%

[X,Y] = meshgrid(t_span_dep,t_span_arr);
close all

% figure(1)
% surf(X,Y,DV);
% colorbar;

% figure(1)
% grid on
% contour(X,Y,DV',100,'LineWidth',1.5)
% colorbar
% caxis([5 10])
% xticks([t_span_dep(1) t_span_dep(360) t_span_dep(744) t_span_dep(1104) t_span_dep(end)])
% xticklabels({'2003 04 01 12:00','2003 05 01 12:00','2003 06 01 12:00','2003 07 01 12:00','2003 08 01 12:00'})
% set(gca,'XTickLabelRotation',45)
% 
% yticks([t_span_arr(1) t_span_arr(360) t_span_arr(732) t_span_arr(1092) t_span_arr(1464) t_span_arr(1836) t_span_arr(end)])
% yticklabels({'2003 09 01 12:00','2003 10 01 12:00','2003 11 01 00','2003 Dec 01','2004 Jan 01','2004 Feb 01','2004 Mar 01'})
% set(gca,'YTickLabelRotation',45)


%%%%%% Plot the transfer arc %%%%%%



TOF = (t_span_arr(t_arr_index) - t_span_dep(t_dep_index)) *24*3600; % [s]
orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
Nrev = 0; %Zero-revolution case
Ncase = 0; %Not used for the zero-revolution case
optionsLMR = 0; %No display
RF_T = RF_T'; % transpose because input of Lambert

[a_T,P,E,ERROR,VI_T,VF,TPAR,THETA] = lambertMR(RI_T,RF_T,TOF,MU,orbitType,Nrev,Ncase,optionsLMR);



%%%%%% Propagation of Earth's orbit %%%%%%

mu_E = astroConstants(13);
RI_E = DEP_eph_E(t_dep_index,1:3); % Earth's position vector at departure
VI_E = DEP_eph_E(t_dep_index,4:6); % Earth's velocity vector at departure

y0_E = [ RI_E'; VI_E']; % transpose because input of ode
T_E = 2*pi*sqrt(a_E^3/ksun);
tspan_E = linspace(0,T_E,1000);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ ~,Y_E ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_E, y0_E, options ); % Y is the position vector of Earth over its orbit


%%%%%%% Propagation of Uranus's orbit %%%%%%

mu_U = astroConstants(17);
RI_U = ARR_eph_U(t_arr_index,1:3); % Jupiter's position vector at departure
VI_U = ARR_eph_U(t_arr_index,4:6); % Jupiter's position vector at departure

y0_U = [ RI_U'; VI_U']; % transpose because input of ode
T_U = 2*pi*sqrt(a_U'^3/ksun);
tspan_U = linspace(0,T_U,1000);

% Perform the integration
[ ~,Y_U ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_U, y0_U, options ); % Y is the position vector of Uranus over its orbit



% Perform the integration of the arc
y0_T = [RI_T; VI_T'];
T_T = TOF;
tspan_T = linspace(0,T_T,1000);
[ ~,Y_T ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_T, y0_T, options ); 

figure 
plotCelBody(0,20e6)
hold on
plot3( Y_E(:,1),Y_E(:,2), Y_E(:,3),'LineWidth',1.5); % Plot of initial orbit (Earth)
plot3(Y_U(:,1),Y_U(:,2),Y_U(:,3),'LineWidth',1.5); % Plot of final orbit (Mars)

plotCelBody(3,8e6,[RI_E(1),RI_E(2),RI_E(3)]) % plot of the Earth (not scaled)
plotCelBody(7,18e6,[RI_U(1),RI_U(2),RI_U(3)]) % plot of Jupiter (not scaled)
% plot3(RI_E(1),RI_E(2),RI_E(3),'ob','MarkerFaceColor','b','MarkerSize',10); % Plot of P1
% plot3(RI_M(1),RI_M(2),RI_M(3),'or','MarkerFaceColor','r','MarkerSize',10) % Plot of P2
axis equal
grid on
comet3(Y_T(:,1),Y_T(:,2),Y_T(:,3)) % Animation of transfer orbit
legend('Sun','Earth''s orbit','Uranus''s orbit','Earth''s position at departure','Uranus''s position at arrival','Transfer arc')
title('Transfer trajectory from Earth to Uranus')
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            %Case 7: Neptune
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'Neptune'

        %%%%%%%%%%%%%%%%%%%%%% Departure time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_dep_in = date2mjd2000([2025 01 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2026 10 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_dep = linspace(t_dep_in,t_dep_fin,1000);% time vector for departure window (it gives me the step every 2 hours)


%%%%%%%%%%%%%%%%%%%%%%% Arrival time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_arr_in = date2mjd2000([2036 01 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2055 06 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_arr = linspace(t_arr_in,t_arr_fin,1000); % time vector for departure window (it gives me the step every 2 hours)



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



%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Mars's ephemerides %%%%%%%%%%%%%%%%%%%%%%

ARR_eph_Nep = [];
mu = ksun;
for k2 = 1:length(t_span_arr)

[kep_arr_Nep,ksun] = uplanet(t_span_arr(k2), 8); % ephemerides of Neptune at arrival window

a_Nep = kep_arr_Nep(1);
e = kep_arr_Nep(2);
i = kep_arr_Nep(3);
Om =  kep_arr_Nep(4);
om =  kep_arr_Nep(5);
f0 =  kep_arr_Nep(6);

[r_arr_Nep,v_arr_Nep] = kep2car_mod(a_Nep,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_arr_Nep = [r_arr_Nep' v_arr_Nep']; % Matrix containing Neptune's state for arrival window


ARR_eph_Nep = [ARR_eph_Nep; s_arr_Nep]; % matrix containing all Mars's states for arrival time window
end

%%
%%%%%%%%%%%%%%%%%%%%% COMPUTING OF LAMBERT ARCS AND DELTA V %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Lambert arcs between departure and arrival windows %%%%%%%%%%%%
    
MU = ksun;
DEP_eph = DEP_eph_E;
ARR_eph = ARR_eph_Nep;
V_inf = 12.5;

[DV_constrained,DV_unconstrained,min_V,date_dep,date_arr,t_dep_index,t_arr_index,DV1,DV2,VI_T,RI_T,RF_T] = Transfer_arc(t_span_dep,t_span_arr,MU,DEP_eph,ARR_eph,V_inf);


%%%%%%%%%%%% Plot of the "Pork-chop plot" %%%%%%%%%%%%%%%%

[X,Y] = meshgrid(t_span_dep,t_span_arr);
close all

% figure(1)
% surf(X,Y,DV);
% colorbar;

% figure(1)
% grid on
% contour(X,Y,DV',100,'LineWidth',1.5)
% colorbar
% caxis([5 10])
% xticks([t_span_dep(1) t_span_dep(360) t_span_dep(744) t_span_dep(1104) t_span_dep(end)])
% xticklabels({'2003 04 01 12:00','2003 05 01 12:00','2003 06 01 12:00','2003 07 01 12:00','2003 08 01 12:00'})
% set(gca,'XTickLabelRotation',45)
% 
% yticks([t_span_arr(1) t_span_arr(360) t_span_arr(732) t_span_arr(1092) t_span_arr(1464) t_span_arr(1836) t_span_arr(end)])
% yticklabels({'2003 09 01 12:00','2003 10 01 12:00','2003 11 01 00','2003 Dec 01','2004 Jan 01','2004 Feb 01','2004 Mar 01'})
% set(gca,'YTickLabelRotation',45)


%%%%%% Plot the transfer arc %%%%%%



TOF = (t_span_arr(t_arr_index) - t_span_dep(t_dep_index)) *24*3600; % [s]
orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
Nrev = 0; %Zero-revolution case
Ncase = 0; %Not used for the zero-revolution case
optionsLMR = 0; %No display
RF_T = RF_T'; % transpose because input of Lambert

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





%%%%%%% Propagation of Neptune's orbit %%%%%%

mu_Nep = astroConstants(18);
RI_Nep = ARR_eph_Nep(t_arr_index,1:3); % Neptune's position vector at departure
VI_Nep = ARR_eph_Nep(t_arr_index,4:6); % Mars's position vector at departure

y0_Nep = [ RI_Nep'; VI_Nep']; % transpose because input of ode
T_Nep = 2*pi*sqrt(a_Nep^3/ksun);
tspan_Nep = linspace(0,T_Nep,1000);

% Perform the integration
[ ~,Y_Nep ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_Nep, y0_Nep, options ); % Y is the position vector of Neptune over its orbit



%%%%%%%%%%%%% Propagation of transfer arc %%%%%%%%%%%%%%%


y0_T = [ RI_T; VI_T']; % transpose because input of ode
T_T = TOF;
tspan_T = linspace(0,T_T,1000);
% Perform the integration
[ ~,Y_T ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_T, y0_T, options ); % Y is the position vector of the transfer arc

%% %%%%%% Plot the transfer arc %%%%%%

figure(2) 
plotCelBody(0,20e6)
hold on
plot3( Y_E(:,1),Y_E(:,2), Y_E(:,3),'LineWidth',1.5); % Plot of initial orbit (Earth)
plot3(Y_Nep(:,1),Y_Nep(:,2),Y_Nep(:,3),'LineWidth',1.5); % Plot of final orbit (Mars)

plotCelBody(3,8e6,[RI_E(1),RI_E(2),RI_E(3)]) % plot of the Earth (not scaled)
plotCelBody(8,4e6,[RI_Nep(1),RI_Nep(2),RI_Nep(3)]) % plot of Neptune (not scaled)
% plot3(RI_E(1),RI_E(2),RI_E(3),'ob','MarkerFaceColor','b','MarkerSize',10); % Plot of P1
% plot3(RI_M(1),RI_M(2),RI_M(3),'or','MarkerFaceColor','r','MarkerSize',10) % Plot of P2
axis equal
grid on
comet3(Y_T(:,1),Y_T(:,2),Y_T(:,3)) % Animation of transfer orbit
legend('Sun','Earth''s orbit','Neptune''s orbit','Earth''s position at departure','Neptune''s position at arrival','Transfer arc')
title('Transfer trajectory from Earth to Neptune')
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')


end