%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ex 1b: Design a flyby around Earth for fixed location of
% the incoming asymptote and different impact parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\timeConversion\time\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\plotCelBody\');


% Data
mu_E = astroConstants(13); % Earth's gravitational constant
mu_Sun = astroConstants(4); % Earth's gravitational constant
v_inf_minus = [15.1; 0; 0]; % [Km/s]
r_E = [astroConstants(2); 0; 0]; % [AU], position vector of Earth wrt Sun
V_planet = [0; sqrt(mu_Sun/astroConstants(2)); 0]; % assuming circular orbit

u_behind = [0; 0; 1]; % I choose this location for the incoming asymptote
u = u_behind;
%%% Compute the 2D hyperbola for different impact parameters %%%%
Impact_p = [9200, 10200, 11200, 12200, 13200]; % [Km]

tspan_1 = linspace(0,-3000);
tspan_2 = linspace(0,3000);

T_E = 365*24*3600; % Earth's period in seconds
tspan_1h = linspace(0,T_E/2);
tspan_2h = linspace(0,-T_E/2);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

Y1_matrix = [];
Y2_matrix = [];
y0_matrix = [];

Y1h_matrix = [];
Y2h_matrix = [];
for i = 1:5

    % Planetocentric frame
    a_hyp = -mu_E/(norm(v_inf_minus)^2); % [Km], semi-major axis
    delta(i) = 2*atan(-a_hyp/Impact_p(i)); % Turn angle
    e_hyp = 1/sin(delta(i)/2); % Eccentricity
    r_p = a_hyp*(1-e_hyp); % [Km], Pericenter radius
    v_p = sqrt(norm(v_inf_minus)^2 + 2*mu_E/r_p); % Pericenter velocity in magnitude
    v_plus = Rodrigues(v_inf_minus,delta(i),u); % v_inf_minus rotated
    deltaV_direction = (v_plus - v_inf_minus)/norm(v_plus - v_inf_minus); % direction of deltaV
    RI = -r_p *  deltaV_direction; % perigee of the hyperbola (vector)
    r_p_direction = RI/(norm(RI)); % Perigee direction
    v_p_direction = cross(u,r_p_direction); % Direction of perigee velocity
    VI = v_p * v_p_direction;
    y0 = [ RI; VI]; % initial condition of hyperbola

   
    % Perform the integration
    [ ~,Y_1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_1, y0, options ); % Y is the position vector of the transfer arc
    [ ~,Y_2 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_2, y0, options ); % Y is the position vector of the transfer arc

    Y1_matrix = [Y1_matrix Y_1]; % matrix containing 6x5 columns (5 = number of impact parameters)
    Y2_matrix = [Y2_matrix Y_2];
    y0_matrix = [y0_matrix y0];

    % Heliocentric frame
    RI_h = [r_E];
    VI_h1 = V_planet + v_plus;
    VI_h2 = V_planet + v_inf_minus;
    y0_h1 = [ RI_h; VI_h1]; % initial condition of hyperbola
    y0_h2 = [ RI_h; VI_h2]; % initial condition of hyperbola
    
    % Perform the integration
    [ ~,Y_1h ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_1h, y0_h1, options ); % Y is the position vector of the transfer arc
    [ ~,Y_2h ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_2h, y0_h2, options ); % Y is the position vector of the transfer arc


    Y1h_matrix = [Y1h_matrix Y_1h];
    Y2h_matrix = [Y2h_matrix Y_2h];
end


%%%%%%% PLOT in planetocentric frame %%%%%%%%%%%%

figure(1)
plotCelBody(3)
hold on
grid on
plot3(Y1_matrix(:,1),Y1_matrix(:,2),Y1_matrix(:,3),'b','LineWidth',1.5); % Plot of the hyperbola 1
plot3(Y2_matrix(:,1),Y2_matrix(:,2),Y2_matrix(:,3),'b','LineWidth',1.5); % Plot of the hyperbola 1
plot3(y0_matrix(1,1),y0_matrix(2,1),y0_matrix(3,1),'ob','MarkerFaceColor','b','MarkerSize',5); % Pericenter point

plot3(Y1_matrix(:,7),Y1_matrix(:,8),Y1_matrix(:,9),'r','LineWidth',1.5); % Plot of the hyperbola 2
plot3(Y2_matrix(:,7),Y2_matrix(:,8),Y2_matrix(:,9),'r','LineWidth',1.5); % Plot of the hyperbola 2
plot3(y0_matrix(1,2),y0_matrix(2,2),y0_matrix(3,2),'or','MarkerFaceColor','r','MarkerSize',5); % Pericenter point

plot3(Y1_matrix(:,13),Y1_matrix(:,14),Y1_matrix(:,15),'y','LineWidth',1.5); % Plot of the hyperbola 3
plot3(Y2_matrix(:,13),Y2_matrix(:,14),Y2_matrix(:,15),'y','LineWidth',1.5); % Plot of the hyperbola 3
plot3(y0_matrix(1,3),y0_matrix(2,3),y0_matrix(3,3),'oy','MarkerFaceColor','y','MarkerSize',5); % Pericenter point

plot3(Y1_matrix(:,19),Y1_matrix(:,20),Y1_matrix(:,21),'m','LineWidth',1.5); % Plot of the hyperbola 4
plot3(Y2_matrix(:,19),Y2_matrix(:,20),Y2_matrix(:,21),'m','LineWidth',1.5); % Plot of the hyperbola 4
plot3(y0_matrix(1,4),y0_matrix(2,4),y0_matrix(3,4),'om','MarkerFaceColor','m','MarkerSize',5); % Pericenter point

plot3(Y1_matrix(:,25),Y1_matrix(:,26),Y1_matrix(:,27),'g','LineWidth',1.5); % Plot of the hyperbola 4
plot3(Y2_matrix(:,25),Y2_matrix(:,26),Y2_matrix(:,27),'g','LineWidth',1.5); % Plot of the hyperbola 4
plot3(y0_matrix(1,5),y0_matrix(2,5),y0_matrix(3,5),'og','MarkerFaceColor','g','MarkerSize',5); % Pericenter point

xlabel('x [R_E]')
ylabel('y [R_E]')
zlabel('z [R_E]')
title('Planetocentric trajectories in Earth-centred frame parallel to HECI, for flybys in front of the planet with same entry velocity and different impact parameters')
legend('','\Delta = 9200 [km]','','','\Delta = 10200 [km]','','','\Delta = 11200 [km]','','','\Delta = 12200 [km]','','','\Delta = 13200 [km]','','Pericenter')


%%%%%%% PLOT in heliocentric frame %%%%%%%%%%%%

figure(2)
plotCelBody(0,695000*20)
hold on
grid on

plot3(Y1h_matrix(:,1),Y1h_matrix(:,2),Y1h_matrix(:,3),'b','LineWidth',1.5); % Plot of the hyperbola 1
plot3(Y2h_matrix(:,1),Y2h_matrix(:,2),Y2h_matrix(:,3),'b','LineWidth',1.5); % Plot of the hyperbola 1

plot3(Y1h_matrix(:,7),Y1h_matrix(:,8),Y1h_matrix(:,9),'r','LineWidth',1.5); % Plot of the hyperbola 2
plot3(Y2h_matrix(:,7),Y2h_matrix(:,8),Y2h_matrix(:,9),'r','LineWidth',1.5); % Plot of the hyperbola 2

plot3(Y1h_matrix(:,13),Y1h_matrix(:,14),Y1h_matrix(:,15),'y','LineWidth',1.5); % Plot of the hyperbola 3
plot3(Y2h_matrix(:,13),Y2h_matrix(:,14),Y2h_matrix(:,15),'y','LineWidth',1.5); % Plot of the hyperbola 3


plot3(Y1h_matrix(:,19),Y1h_matrix(:,20),Y1h_matrix(:,21),'m','LineWidth',1.5); % Plot of the hyperbola 4
plot3(Y2h_matrix(:,19),Y2h_matrix(:,20),Y2h_matrix(:,21),'m','LineWidth',1.5); % Plot of the hyperbola 4

plot3(Y1h_matrix(:,25),Y1h_matrix(:,26),Y1h_matrix(:,27),'g','LineWidth',1.5); % Plot of the hyperbola 4
plot3(Y2h_matrix(:,25),Y2h_matrix(:,26),Y2h_matrix(:,27),'g','LineWidth',1.5); % Plot of the hyperbola 4


xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')
title('Helioocentric trajectories in HECI frame, for flybys in front of the planet with same entry velocity and different impact parameters')
legend('','\Delta = 9200 [km]','','\Delta = 10200 [km]','','\Delta = 11200 [km]','','\Delta = 12200 [km]','','\Delta = 13200 [km]','')
