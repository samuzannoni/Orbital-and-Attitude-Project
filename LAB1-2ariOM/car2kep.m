function [a, e, i, OM, om, th] = car2kep(r,v,mu)

% car2kep.m - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
% [a, e, i, OM, om, th] = car2kep(r, v, mu)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in
% radians.
%
% INPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% i [1x1] Inclination [rad]
% OM [1x1] RAAN [rad]
% om [1x1] Pericentre anomaly [rad]
% th [1x1] True anomaly [rad]

R=norm(r);
V=norm(v);
I=[1 0 0]';
J=[0 1 0]';
K=[0 0 1]';

a=1/(2 / R - V^2/mu);  

h=cross(r,v);
H=norm(h);

e=(cross(v,h) ./ mu) - (r./R);
E=norm(e);

i=acos(h(3,1)/H);

n=cross(K,h);
N=norm(n);


if n(2,1)>= 0
    OM=acos(dot(n,I)/N);
else 
    OM=2*pi -acos(dot(n,I)/N);
end

if dot(e,K) >= 0
    om=acos(dot(n,e)/(N*E));
else
    om=2*pi - acos(dot(n,e)/(N*E));
end

if dot(r,v)>=0
    th=acos(dot(r,e)/(R*E));
else
    th=2*pi-acos(dot(r,e)/(R*E));
end
e=E;
end
