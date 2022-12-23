function [rp_hyp,a_hyp,e_hyp,deltadegree,deltav_pnorm,delta] = hyperbola(vinfminus,ip,mu)
%%% INPUT
% ip: impact parameter
% vinfminus: excess velocity(minus)
%%% OUTPUT
% rp_hyp: radius of pericentre
% delta: turning angle

vinfnorm=norm(vinfminus);
a_hyp = -mu/(vinfnorm.^2);
delta = 2*atan(-a_hyp./ip);
e_hyp = 1./sin(delta./2);
rp_hyp= a_hyp.*(1-e_hyp);
deltadegree = rad2deg(delta);

%beta = (pi-delta)./2;

deltav_pnorm= 2*vinfnorm*sin(delta./2);

end