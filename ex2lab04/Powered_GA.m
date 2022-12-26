function [rp] = PGA(turn_angle,v_inf_minus,v_inf_plus,mu,rp0)
vinfminus= norm(v_inf_minus);
vinfplus= norm(v_inf_plus);

turn_angle_solve = @(rp) asin(1./(1+ rp*vinfminus^2/mu))+asin(1./(1+ rp*vinfplus^2/mu))-turn_angle;

options = optimoptions('fsolve','TolFun',1e-14,'Display','off');

rp = fsolve(turn_angle_solve, rp0, options);
end
