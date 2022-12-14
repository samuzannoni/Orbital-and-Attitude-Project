clc
clear all
close all

addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\timeConversion\time\');
addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\plotCelBody\');

%%%%% Data %%%%%
mu_E = astroConstants(13); % Earth's gravitational constant
mu_Sun = astroConstants(4); % Earth's gravitational constant
v_inf_minus = [15.1; 0; 0]; % [Km/s]
Impact_p = 9200; % [Km]
r_E = astroConstants(2); 
%%% Compute the hyperbola %%%%
a_hyp = -mu_E/(norm(v_inf_minus)^2); % [Km], semi-major axis
delta = 2*atan(-a_hyp/Impact_p); % Turn angle
e_hyp = 1/sin(delta/2); % Eccentricity
r_p = a_hyp*(1-e_hyp); % [Km], Pericenter radius 

V_planet = [0; sqrt(mu_Sun/astroConstants(2)); 0]; % assuming circular orbit

%%%% Computing v_inf_plus for 3 different cases: in front, behind and under the planet %%%%%

u_front = [0; 0; -1];
u_behind = [0; 0; 1];
u_under = [0; -1; 0];

u = [u_front u_behind u_under];

% V_inf_plus = [];
% 
% for i = 1:size(u)
% 
%     V_rot = v_inf_minus * cos(delta) + cross(u(:,i),v_inf_minus) .* sin(delta) + (u(:,1)*(dot(u(:,i),v_inf_minus))) * (1 - cos(delta));
%    
%     V_inf_plus = [V_inf_plus V_rot]; % Matrix containing v_inf_plus (each column corresponds to a different case)
% end
% 
% %%%%% Compute Delta V %%%%%%
% 
% V_E = sqrt(mu_Sun/r_E);

%%%  Fly-by in front of the planet %%%%

u = u_front;
v_plus_f = Rodrigues(v_inf_minus,delta,u); % excess velocity in planetocentric frame
V_minus_f = V_planet + v_inf_minus; % Incoming velocity wrt the Sun
V_plus_f = V_planet + v_plus_f; % Excess velocity wrt the Sun
%%%  Fly-by behind of the planet %%%%

u = u_behind;
v_plus_b = Rodrigues(v_inf_minus,delta,u); % excess velocity in planetocentric frame
V_minus_b = V_planet + v_inf_minus; % Incoming velocity wrt the Sun
V_plus_b = V_planet + v_plus_b; % Excess velocity wrt the Sun

%%%  Fly-by behind of the planet %%%%

u = u_under;
v_plus_u = Rodrigues(v_inf_minus,delta,u); % excess velocity in planetocentric frame
V_minus_u = V_planet + v_inf_minus; % Incoming velocity wrt the Sun
V_plus_u = V_planet + v_plus_u; % Excess velocity wrt the Sun


%%%%%%% Propagation of hyperbolas %%%%%%%%%%%%

Flyby = input('Choose Flyby :');

 tspan_1 = linspace(0,-3000);
 tspan_2 = linspace(0,3000);

switch Flyby

case 'in front'

    close all

    u = u_front;
    %%%% Planetocentric frame %%%%

    deltaV_direction = (v_plus_f - v_inf_minus)/norm(v_plus_f - v_inf_minus); % direction of deltaV
    RI_front = -r_p *  deltaV_direction; % perigee of the hyperbola (vector)
    r_p_direction = RI_front/(norm(RI_front)); % Perigee direction
    v_p = sqrt(norm(v_inf_minus)^2 + 2*mu_E/r_p); % Pericenter velocity in magnitude
    v_p_direction = cross(u,r_p_direction); % Direction of perigee velocity
    VI_front = v_p * v_p_direction;

    y0_front = [ RI_front; VI_front]; % initial condition of hyperbola


    options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
    % Perform the integration
    [ ~,Y_front1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_1, y0_front, options ); % Y is the position vector of the transfer arc
    [ ~,Y_front2 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_2, y0_front, options ); % Y is the position vector of the transfer arc

    %%%%% PLOT  in planetocentric %%%%%%%%

    figure(1)
    plotCelBody(3); % Earth plot
    hold on
    plot3(Y_front1(:,1),Y_front1(:,2),Y_front1(:,3),'b','LineWidth',1.5); % Plot of the hyperbola
    plot3(Y_front2(:,1),Y_front2(:,2),Y_front2(:,3),'b','LineWidth',1.5); % Plot of the hyperbola
    plot3(y0_front(1),y0_front(2),y0_front(3),'ob','MarkerFaceColor','b','MarkerSize',5); % Pericenter point
    xlabel('x [R_E]');
    ylabel('y [R_E]');
    zlabel('z [R_E]');
    title('Flyby in front of the planet')
    grid on

    %%%% PLOT in heliocentric frame %%%%%

    RI_fronth = [astroConstants(2); 0; 0]; % radius of the Earth's orbit (circular)
    VI_fronth1 = V_plus_f;
    y0_fronth1 = [ RI_fronth; VI_fronth1];

    VI_fronth2 = V_minus_f;
    y0_fronth2 = [ RI_fronth; VI_fronth2];

    % Perform the integration
    T_E = 365*24*3600; % Earth's period in seconds
    tspan_1h = linspace(0,T_E/2);
    tspan_2h = linspace(0,-T_E/2);

    [ ~,Y_fronth1 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_1h, y0_fronth1, options ); % Y is the position vector of the transfer arc
    [ ~,Y_fronth2 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_2h, y0_fronth2, options ); % Y is the position vector of the transfer arc

    figure(2)
    plotCelBody(0,695000*10');
    hold on
    
    plot3(Y_fronth1(:,1),Y_fronth1(:,2),Y_fronth1(:,3),'r','LineWidth',1.5); % Plot of the hyperbola
    
    plot3(Y_fronth2(:,1),Y_fronth2(:,2),Y_fronth2(:,3),'b','LineWidth',1.5); % Plot of the hyperbola


    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    title('Trajectory in heliocentric frame (HECI) frame (leading-side flyby)')
    grid on
    legend('','After flyby','before flyby')


case 'behind'

    close all

    u = u_behind;
    
    deltaV_direction = (v_plus_b - v_inf_minus)/norm(v_plus_b - v_inf_minus); % direction of deltaV
    RI_behind = -r_p *  deltaV_direction; % perigee of the hyperbola (vector)
    r_p_direction = RI_behind/(norm(RI_behind)); % Perigee direction
    v_p = sqrt(norm(v_inf_minus)^2 + 2*mu_E/r_p); % Pericenter velocity in magnitude
    v_p_direction = cross(u,r_p_direction); % Direction of perigee velocity
    VI_behind = v_p * v_p_direction;

    y0_behind = [ RI_behind; VI_behind]; % initial condition of hyperbola
   

    options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
    % Perform the integration
    [ ~,Y_behind1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_1, y0_behind, options ); % Y is the position vector of the transfer arc
    [ ~,Y_behind2 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_2, y0_behind, options ); % Y is the position vector of the transfer arc

    %%%% PLOT: Planetocentric frame %%%%%%%%%%%%

    plotCelBody(3); % Earth plot
    hold on
    plot3(Y_behind1(:,1),Y_behind1(:,2),Y_behind1(:,3),'g','LineWidth',1.5); % Plot of the hyperbola
    plot3(Y_behind2(:,1),Y_behind2(:,2),Y_behind2(:,3),'g','LineWidth',1.5); % Plot of the hyperbola
    plot3(y0_behind(1),y0_behind(2),y0_behind(3),'og','MarkerFaceColor','g','MarkerSize',5); % Pericenter point
    xlabel('x [R_E]');
    ylabel('y [R_E]');
    zlabel('z [R_E]');
    title('Flyby behind the planet')
    grid on


    %%%%% PLOT in heliocentric frame %%%%%

    RI_behind_h = [astroConstants(2); 0; 0]; % radius of the Earth's orbit (circular)
    VI_behind_h1 = V_plus_b;
    y0_behind_h1 = [ RI_behind_h; VI_behind_h1];

    VI_behind_h2 = V_minus_b;
    y0_behind_h2 = [ RI_behind_h; VI_behind_h2];

    % Perform the integration
    T_E = 365*24*3600; % Earth's period in seconds
    tspan_1h = linspace(0,T_E/2);
    tspan_2h = linspace(0,-T_E/2);

    [ ~,Y_behind_h1 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_1h, y0_behind_h1, options ); % Y is the position vector of the transfer arc
    [ ~,Y_behind_h2 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_2h, y0_behind_h2, options ); % Y is the position vector of the transfer arc

    figure(2)
    plotCelBody(0,695000*10');
    hold on
    
    plot3(Y_behind_h1(:,1),Y_behind_h1(:,2),Y_behind_h1(:,3),'r','LineWidth',1.5); % Plot of the hyperbola
    
    plot3(Y_behind_h2(:,1),Y_behind_h2(:,2),Y_behind_h2(:,3),'b','LineWidth',1.5); % Plot of the hyperbola


    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    title('Trajectory in heliocentric frame (HECI) frame (trailing-side flyby)')
    grid on
    legend('','After flyby','before flyby')

case 'under'

    close all

        u = u_under;

    %%%% Planetocentric frame %%%%
    
    deltaV_direction = (v_plus_u - v_inf_minus)/norm(v_plus_u - v_inf_minus); % direction of deltaV
    RI_under = -r_p *  deltaV_direction; % perigee of the hyperbola (vector)
    r_p_direction = RI_under/(norm(RI_under)); % Perigee direction
    v_p = sqrt(norm(v_inf_minus)^2 + 2*mu_E/r_p); % Pericenter velocity in magnitude
    v_p_direction = cross(u,r_p_direction); % Direction of perigee velocity
    VI_under = v_p * v_p_direction;

    y0_under = [ RI_under; VI_under]; % initial condition of hyperbola
    

    options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
    % Perform the integration
    [ ~,Y_under1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_1, y0_under, options ); % Y is the position vector of the transfer arc
    [ ~,Y_under2 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_2, y0_under, options ); % Y is the position vector of the transfer arc

    %%%% PLOT %%%%%%%%

    plotCelBody(3); % Earth plot
    hold on
    plot3(Y_under1(:,1),Y_under1(:,2),Y_under1(:,3),'r','LineWidth',1.5); % Plot of the hyperbola
    plot3(Y_under2(:,1),Y_under2(:,2),Y_under2(:,3),'r','LineWidth',1.5); % Plot of the hyperbola
    plot3(y0_under(1),y0_under(2),y0_under(3),'or','MarkerFaceColor','r','MarkerSize',5);
    xlabel('x [R_E]');
    ylabel('y [R_E]');
    zlabel('z [R_E]');
    title('Flyby under the planet')
    grid on



     %%%%% PLOT in heliocentric frame %%%%%

    RI_under_h = [astroConstants(2); 0; 0]; % radius of the Earth's orbit (circular)
    VI_under_h1 = V_plus_u;
    y0_under_h1 = [ RI_under_h; VI_under_h1];

    VI_under_h2 = V_minus_b;
    y0_under_h2 = [ RI_under_h; VI_under_h2];

    % Perform the integration
    T_E = 365*24*3600; % Earth's period in seconds
    tspan_1h = linspace(0,T_E/2);
    tspan_2h = linspace(0,-T_E/2);

    [ ~,Y_under_h1 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_1h, y0_under_h1, options ); % Y is the position vector of the transfer arc
    [ ~,Y_under_h2 ] = ode113( @(t,y) ode_2bp(t,y,mu_Sun), tspan_2h, y0_under_h2, options ); % Y is the position vector of the transfer arc

    figure(2)
    plotCelBody(0,695000*10');
    hold on
    
    plot3(Y_under_h1(:,1),Y_under_h1(:,2),Y_under_h1(:,3),'r','LineWidth',1.5); % Plot of the hyperbola
    
    plot3(Y_under_h2(:,1),Y_under_h2(:,2),Y_under_h2(:,3),'b','LineWidth',1.5); % Plot of the hyperbola


    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    title('Trajectory in heliocentric frame (HECI) frame (under Earth flyby)')
    grid on
    legend('','After flyby','before flyby')

    case 'All flyby'

        close all

        u = u_front;
        %%%% Planetocentric frame %%%%

        deltaV_direction = (v_plus_f - v_inf_minus)/norm(v_plus_f - v_inf_minus); % direction of deltaV
        RI_front = -r_p *  deltaV_direction; % perigee of the hyperbola (vector)
        r_p_direction = RI_front/(norm(RI_front)); % Perigee direction
        v_p = sqrt(norm(v_inf_minus)^2 + 2*mu_E/r_p); % Pericenter velocity in magnitude
        v_p_direction = cross(u,r_p_direction); % Direction of perigee velocity
        VI_front = v_p * v_p_direction;

        y0_front = [ RI_front; VI_front]; % initial condition of hyperbola


        options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
        % Perform the integration
        [ ~,Y_front1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_1, y0_front, options ); % Y is the position vector of the transfer arc
        [ ~,Y_front2 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_2, y0_front, options ); % Y is the position vector of the transfer arc

        plotCelBody(3); % Earth plot
        hold on
        plot3(Y_front1(:,1),Y_front1(:,2),Y_front1(:,3),'b','LineWidth',1.5); % Plot of the hyperbola
        plot3(Y_front2(:,1),Y_front2(:,2),Y_front2(:,3),'b','LineWidth',1.5); % Plot of the hyperbola
        plot3(y0_front(1),y0_front(2),y0_front(3),'ob','MarkerFaceColor','b','MarkerSize',5);




        u = u_behind;


        %%%% Planetocentric frame %%%%

        deltaV_direction = (v_plus_b - v_inf_minus)/norm(v_plus_b - v_inf_minus); % direction of deltaV
        RI_behind = -r_p *  deltaV_direction; % perigee of the hyperbola (vector)
        r_p_direction = RI_behind/(norm(RI_behind)); % Perigee direction
        v_p = sqrt(norm(v_inf_minus)^2 + 2*mu_E/r_p); % Pericenter velocity in magnitude
        v_p_direction = cross(u,r_p_direction); % Direction of perigee velocity
        VI_behind = v_p * v_p_direction;

        y0_behind = [ RI_behind; VI_behind]; % initial condition of hyperbola


        options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
        % Perform the integration
        [ ~,Y_behind1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_1, y0_behind, options ); % Y is the position vector of the transfer arc
        [ ~,Y_behind2 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_2, y0_behind, options ); % Y is the position vector of the transfer arc




        plot3(Y_behind1(:,1),Y_behind1(:,2),Y_behind1(:,3),'g','LineWidth',1.5); % Plot of the hyperbola
        plot3(Y_behind2(:,1),Y_behind2(:,2),Y_behind2(:,3),'g','LineWidth',1.5); % Plot of the hyperbola
        plot3(y0_behind(1),y0_behind(2),y0_behind(3),'og','MarkerFaceColor','g','MarkerSize',5); % Pericenter point



        u = u_under;

        %%%% Planetocentric frame %%%%

        deltaV_direction = (v_plus_u - v_inf_minus)/norm(v_plus_u - v_inf_minus); % direction of deltaV
        RI_under = -r_p *  deltaV_direction; % perigee of the hyperbola (vector)
        r_p_direction = RI_under/(norm(RI_under)); % Perigee direction
        v_p = sqrt(norm(v_inf_minus)^2 + 2*mu_E/r_p); % Pericenter velocity in magnitude
        v_p_direction = cross(u,r_p_direction); % Direction of perigee velocity
        VI_under = v_p * v_p_direction;

        y0_under = [ RI_under; VI_under]; % initial condition of hyperbola


        options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
        % Perform the integration
        [ ~,Y_under1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_1, y0_under, options ); % Y is the position vector of the transfer arc
        [ ~,Y_under2 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_2, y0_under, options ); % Y is the position vector of the transfer arc

        %%%% PLOT %%%%%%%%

        plotCelBody(3); % Earth plot
        plot3(Y_under1(:,1),Y_under1(:,2),Y_under1(:,3),'r','LineWidth',1.5); % Plot of the hyperbola
        plot3(Y_under2(:,1),Y_under2(:,2),Y_under2(:,3),'r','LineWidth',1.5); % Plot of the hyperbola
        plot3(y0_under(1),y0_under(2),y0_under(3),'or','MarkerFaceColor','r','MarkerSize',5);


        xlabel('x [R_E]');
        ylabel('y [R_E]');
        zlabel('z [R_E]');
        title('Flyby in Earth-centred frame parallel to the heliocentric inertial (HECI) frame')
   
        legend('','Flyby in front of the planet','','','Flyby behind the planet','','','Flyby under the planet','','Perigee');

 end

