function dy = ode_2bp_j2( ~, y, mu_E,R_Earth,J2 )
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3);
v = y(4:6);

r_norm = norm(r);

a_j2 = 3/2 * (J2*mu_E*R_Earth^2)/r_norm^4 * [r(1)/r_norm * (5 * r(3)^2/r_norm^2 - 1), r(2)/r_norm * (5 * r(3)^2/r_norm^2 - 1), r(3)/r_norm * (5 * r(3)^2/r_norm^2 - 3) ]; 

% Set the derivatives of the state
dy = [ v(1), v(2), v(3), (-(mu_E/r_norm^3)*r(1) + a_j2(1)), (-(mu_E/r_norm^3)*r(2)+a_j2(2)), (-(mu_E/r_norm^3)*r(3) + a_j2(3)) ]';
end
