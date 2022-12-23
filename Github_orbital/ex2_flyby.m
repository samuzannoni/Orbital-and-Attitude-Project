%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design a powered gravity assist around the Earth given the 
% heliocentric velocities bedore and after the flyby, and
% Earth's position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all



addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\timeConversion\time\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\plotCelBody\');



%%%%%%% DATA %%%%%%%%%%%

mu_E = astroConstants(13); % Earth's Gravitational constant
mu_sun = astroConstants(4); % Sun's Gravitational constant
r_E = [0; -astroConstants(2); 0]; % [AU], position vector of Earth wrt Sun
V_minus = [31.5; 5.2; 0.0]; % [Km/s], heliocentric entry velocity
V_plus = [36.0; 0.0; 0.0]; % [Km/s], heliocentric exit velocity

%%
%%%% 1: Compute the velocities relative to the planet before and after the flyby %%%%%%

V_Earth = [sqrt(mu_sun/astroConstants(2)); 0.0; 0.0]; % Velocity of the Earth
v_inf_minus = V_minus - V_Earth; % entry velocity (of the s/c)
v_inf_plus = V_plus - V_Earth; % exit velocity (of the s/c)

%%
%%%% 2 & 3: Compute the turning angle and rp and check feasibility %%%%%%%%%%%%%%%%%%%%%%%%%

rp0 = astroConstants(23); % Earth's radius
MU = mu_E;

turn_angle = acos(v_inf_plus'*v_inf_minus/(norm(v_inf_plus)*norm(v_inf_minus)));
delta = rad2deg(turn_angle);

[rp] = Powered_GA(v_inf_minus,v_inf_plus,turn_angle,rp0,MU);
delta_minus_half = 2*asin(1/(1 + (rp*norm(v_inf_minus)^2/MU)))/2;
delta_plus_half = 2*asin(1/(1 + (rp*norm(v_inf_plus)^2/MU)))/2;
h_ga = rp - astroConstants(23); % altitude of the hyperbola
h_atm = 100; % [Km], Atmosphere's altitude (Karman line) 

fprintf('The total turn angle is %.6g deg \n',delta)
fprintf('The radius of the pericenter is %.8g [Km]\n',rp)

if h_ga > h_atm

    fprintf('The gravity assist is performed at %.8g [Km] \n',h_ga)

else

    fprintf('The altitude is not enough to perform a flyby')

end

%%%%%% 4: Compite the velocities of the two hyperbolic arcs at percienter and the required Delta V 

vp_entry = sqrt(norm(v_inf_minus)^2 + 2*mu_E/rp); % Pericenter velocity in magnitude of entry hyperbola
vp_exit = sqrt(norm(v_inf_plus)^2 + 2*mu_E/rp); % Pericenter velocity in magnitude of exit hyperbola
Delta_vp = vp_entry - vp_exit;
Delta_vp_norm = norm(vp_exit - vp_entry);
fprintf('The DeltaV_p is %.5g [Km/s] \n',Delta_vp_norm);


%%%%%%% 5: Plot in planetocentric and heliocentric frame %%%%%%%%
u = (cross(v_inf_minus,v_inf_plus)) / norm(cross(v_inf_minus,v_inf_plus));

v_plus = Rodrigues(v_inf_minus,delta,u); % v_inf_minus rotated (considering no change in modulus)
deltaV_direction = (v_plus - v_inf_minus)/norm(v_plus - v_inf_minus); % direction of deltaV (free flyby)
vp_direction = -cross(u,deltaV_direction); % Direction of delta V at pericenter
rp_direction = -cross(u,vp_direction);
% Beta_minus = pi/2 - delta_minus_half; % [Rad], beta angle
% rp_direction = (v_inf_minus * cos(Beta_minus))/(norm(v_inf_minus * cos(Beta_minus))); % direction of the pericenter
Rp_vector = rp * rp_direction; % rp vector, initial condition to ode
% vp_direction = (v_inf_minus * sin(Beta_minus))/ (norm( v_inf_minus * sin(Beta_minus)));
Vp_vector_entry = vp_entry * vp_direction;  % vp entry vector, initial condition to first ode
Vp_vector_exit = vp_exit * vp_direction; % vp exit vector, initial condition to second ode

y0_entry = [ Rp_vector; Vp_vector_entry]; % initial condition of entry hyperbola
y0_exit = [ Rp_vector; Vp_vector_exit]; % initial condition of entry hyperbola

tspan_1 = linspace(0,-7000);
tspan_2 = linspace(0,7000);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ ~,Y_entry ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_1, y0_entry, options ); 
[ ~,Y_exit ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_2, y0_exit, options ); 

 %%%%% PLOT  in planetocentric %%%%%%%%

 figure(1)
 plotCelBody(3); % Earth plot
 hold on
 plot3(Y_entry(:,1),Y_entry(:,2),Y_entry(:,3),'b','LineWidth',1.5); % Plot of the hyperbola
 plot3(Y_exit(:,1),Y_exit(:,2),Y_exit(:,3),'r','LineWidth',1.5); % Plot of the hyperbola
 plot3(y0_entry(1),y0_entry(2),y0_entry(3),'ob','MarkerFaceColor','b','MarkerSize',5); % Pericenter point
 xlabel('x [R_E]');
 ylabel('y [R_E]');
 zlabel('z [R_E]');
 title('Powered gravity assist (planetocentric frame)')
 legend('','Incoming hyperbola','Outcoming hyperbola')
 grid on

 %%%%%% 6: Compute all the parameters of the two hyperbolas %%%%%%%%%%

 e_hyp_entry = 1 + (rp*(norm(v_inf_minus)^2)/mu_E); % Eccentricity of entry hyperbola
 e_hyp_exit = 1 + (rp*(norm(v_inf_plus)^2)/mu_E); % Eccentricity of exit hyperbola
 a_hyp_entry = -mu_E/(norm(v_inf_minus)^2); % [Km], semi-major axis of entry hyperbola
 a_hyp_exit = -mu_E/(norm(v_inf_plus)^2); % [Km], semi-major axis of exit hyperbola
 deltaV_flyby = (v_plus - v_inf_minus); % Delta V due to the flyby;


fprintf('The entry eccentricity is %.6g\n',e_hyp_entry)
fprintf('The exit eccentricity is %.6g \n',e_hyp_exit)
fprintf('The entry semi-major axis is %.8g [Km]\n',a_hyp_entry)
fprintf('The exit semi-major axis is %.8g [Km]\n',a_hyp_exit)
fprintf('The DeltaV given by the flyby is %.5g [Km/s] \n',norm(deltaV_flyby))



%%%%% Compute it in heliocentric frame %%%%%%%%


% Heliocentric frame

RI_h = [r_E];
VI_h1 = V_plus;
VI_h2 = V_minus;
y0_h1 = [ RI_h; VI_h1]; % initial condition of hyperbola
y0_h2 = [ RI_h; VI_h2]; % initial condition of hyperbola

T_E = 365*24*3600; % Earth's period in seconds
tspan_1h = linspace(0,T_E/2);
tspan_2h = linspace(0,-T_E/2);


% Perform the integration
[ ~,Y_1h ] = ode113( @(t,y) ode_2bp(t,y,mu_sun), tspan_1h, y0_h1, options ); % Y is the position vector of the transfer arc
[ ~,Y_2h ] = ode113( @(t,y) ode_2bp(t,y,mu_sun), tspan_2h, y0_h2, options ); % Y is the position vector of the transfer arc


 figure(2)
 plotCelBody(0,695000*10); % Sun plot
 hold on
 plot3(Y_1h(:,1),Y_1h(:,2),Y_1h(:,3),'b','LineWidth',1.5); % Plot of the hyperbola
 plot3(Y_2h(:,1),Y_2h(:,2),Y_2h(:,3),'r','LineWidth',1.5); % Plot of the hyperbola
 plot3(y0_h1(1),y0_h1(2),y0_h1(3),'ob','MarkerFaceColor','b','MarkerSize',5); % Pericenter point
 xlabel('x [AU]');
 ylabel('y [AU]');
 zlabel('z [AU]');
 title('Powered gravity assist (heliocentric frame)')
 legend('','Incoming hyperbola','Outcoming hyperbola')
 grid on

 