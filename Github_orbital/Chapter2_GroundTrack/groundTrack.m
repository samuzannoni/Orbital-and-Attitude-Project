function [alpha,delta,lon,lat] = groundTrack(r0,v0,theta_G0,tspan_orbit,w_E,mu_E,t0)

%% orbit propagation

y0 = [ r0; v0 ];
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan_orbit, y0, options );

%% Conversion from cartesian coordinate to RA and Dec

r = Y(:,1:3);
r_norm = vecnorm(r,2,2);
delta = asin(Y(:,3)./r_norm); % Declination

% if (Y(:,2)/r_norm) > 0
%     alpha = acos (Y(:,1)./ (r_norm .* cos(delta)));
% else 
%      alpha = 2*pi - acos (Y(:,1)./ (r_norm .* cos(delta)));
% end

alpha = atan2(Y(:,2),Y(:,1)) % Right ascension


%% Conversion to longitude and latitude


theta_G = @(t) theta_G0 + w_E * (t - t0);

lon = alpha' - theta_G(tspan_orbit); % longitude
lon = lon';
lat = delta; % latitude

end

