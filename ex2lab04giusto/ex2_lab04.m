%% FLYBY 2
% %DATA 

mu_E = astroConstants(13); % Earth's Gravitational constant
mu_sun = astroConstants(4); % Sun's Gravitational constant
r_E = [0; -astroConstants(2); 0]; % [AU], position vector of Earth wrt Sun
V_minus = [31.5; 5.2; 0.0]; % [Km/s], heliocentric entry velocity
V_plus = [36.0; 0.0; 0.0]; % [Km/s], heliocentric exit velocity

%% 1: Compute the velocities relative to the planet before and after the flyby %%%%%%

V_Earth = [sqrt(mu_sun/astroConstants(2)); 0; 0]; % Velocity of the Earth

v_inf_minus = V_minus - V_Earth; % entry velocity (of the s/c)
v_inf_plus = V_plus - V_Earth; % exit velocity (of the s/c)

vinfminus= norm(v_inf_minus);
vinfplus= norm(v_inf_plus);
%% 2: Compute the turning angle 
turn_angle = acos(dot(v_inf_minus,v_inf_plus)/(vinfminus*vinfplus));
turn_angle_deg = rad2deg(turn_angle);

%% 3: Solve the non-linear system for 𝑟𝑝 and check its validity.

rp0 = astroConstants(23); % Earth's radius

[rp] = PGA(turn_angle,v_inf_minus,v_inf_plus,mu_E,rp0);  % Periapsis radiu of the hyperbolas [km]

h_ga = rp - rp0; % altitude of the hyperbola
h_atm = 100; % [Km], Atmosphere's altitude (Karman line) 

if h_ga > h_atm

    disp('The altitude is enough to perform a flyby')

else

    disp('The altitude is not enough to perform a flyby')

end

%% 4: Compute the velocities of the two hyperbolic arcs at pericentre and the required 𝛥𝑣p
v_p_minus = sqrt(vinfminus^2 + 2*mu_E/rp);
v_p_plus = sqrt(vinfplus^2 + 2*mu_E/rp);

deltav_p_PGA = abs(v_p_plus - v_p_minus);
%%
h_ga = rp - rp0;

delta_minus  =  2*asin(1./(1+ rp*vinfminus^2/mu_E));
delta_plus   =  2*asin(1./(1+ rp*vinfplus^2/mu_E));
e_minus= 1+rp*vinfminus^2/mu_E;
e_plus= 1+rp*vinfplus^2/mu_E;
a_minus = -mu_E/(vinfminus.^2);
a_plus = -mu_E/(vinfplus.^2);

delta_V_flyby= norm(v_inf_plus-v_inf_minus);

%% Plot the trajectories (planetocentric)

u = cross(v_inf_plus,v_inf_minus) ./ norm(cross(v_inf_plus,v_inf_minus));

betaplus =(pi-delta_plus)/2;

rp_direction = v_inf_plus./vinfplus*cos(betaplus)./norm(v_inf_plus./vinfplus*cos(betaplus));
vp_direction = cross(u,rp_direction)/norm(cross(u,rp_direction))  ;% Direction of delta V at pericenter

rpvector=rp.*rp_direction;

vp_minus = v_p_minus.*vp_direction;
vp_plus = v_p_plus.*vp_direction;

y0minus = [rpvector;vp_minus];
y0plus = [rpvector;vp_plus];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

tspan1=linspace(0,-7000,100);
tspan2=linspace(0,7000,100);

[TE, Y_hyp1] = ode113( @(t,y) ode_2bp(t,y, mu_E), tspan1, y0minus, options );
[TE, Y_hyp2] = ode113( @(t,y) ode_2bp(t,y, mu_E), tspan2, y0plus, options );

figure()
plot3( Y_hyp1(:,1), Y_hyp1(:,2), Y_hyp1(:,3),'b','LineWidth',2 );
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
axis equal;
grid on;
hold on;
plot3( Y_hyp2(:,1), Y_hyp2(:,2), Y_hyp2(:,3),'r','LineWidth',2 );


asymptote_minus = @(x)   tan(pi/2 - delta_minus/2) * (x -(rp+abs(a_minus)));
asymptote_plus =  @(x)  -tan(pi/2 - delta_plus/2)  * (x -(rp+abs(a_plus)));

X1 = [min(Y_hyp1(:,1)):(rp+abs(a_minus))];
X2 = [min(Y_hyp2(:,1)):(rp+abs(a_plus))];

plot3(X1, asymptote_minus(X1), zeros(length(X1),1), '--k', 'LineWidth', .5)
plot3(X2, asymptote_plus(X2), zeros(length(X2),1), '--g', 'LineWidth', .5)

legend('incoming hyp', 'outcoming hyp','incoming asymptote','outcoming asymptote','Earth');
Terra3d
%% heliocentric
y0_helio1 = [r_E;V_minus];
y0_helio2 = [r_E;V_plus];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

tspan1 = linspace(0,-365*24*3600/2,100);
tspan2=linspace(0,365*24*3600/2,100);

% Perform the integration
[TE, Y_helio1] = ode113( @(t,y) ode_2bp(t,y, mu_sun), tspan1, y0_helio1, options );
[TM, Y_helio2] = ode113( @(t,y) ode_2bp(t,y, mu_sun), tspan2, y0_helio2, options );

%Plot
figure
plot3( Y_helio1(:,1), Y_helio1(:,2), Y_helio1(:,3),'b' );
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
axis equal;
grid on;
hold on;
plot3( Y_helio2(:,1), Y_helio2(:,2), Y_helio2(:,3),'r' );
legend('before flyby', 'after flyby');


