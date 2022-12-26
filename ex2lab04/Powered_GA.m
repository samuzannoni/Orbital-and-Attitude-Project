[rp] = = Powered_GA(v_inf_minus,v_inf_plus,turn_angle,rp0,mu);
delta_minus  = @(rp) 2*asin(1./(1+ rp.*norm(v_inf_minus)^2/mu));
delta_plus   = @(rp) 2*asin(1./(1+ rp.*norm(v_inf_plus)^2/mu));
r_p_SolveFun = @(rp) (delta_minus(rp) + delta_plus(rp))/2 - delta;
r_p_min = R_planet + h_atm;
options = optimoptions('fsolve','TolFun',1e-14,'Display','off');

rp = fsolve(r_p_SolveFun, r_p_min, options);  % Periapsis radiu of the hyperbolas [km]


% rp0 = astroConstants(23);
% delta = @(rp) asin(1/(1+rp.*norm(v_inf_minus)^2/mu_E))+asin(1/(1+rp.*norm(v_inf_plus)^2/mu_E));
% rp = fzero(delta,rp0);