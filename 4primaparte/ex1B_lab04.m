%% 1.B
mu_E = astroConstants(13);
mu_S = astroConstants(4);
AU = astroConstants(2);
vinfminus=[15.1; 0; 0];
r_E=[1; 0; 0]*AU;

uset= ufront; % choose a location for the incoming asymptote
ip = linspace(9200,13200,5);

rp_hyp_=[];
a_hyp_=[];
e_hyp_=[];
deltav_pnorm_=[];
delta_=[];

rpvector_=[];
v_rot_=[];
for k=1:length(ip)
    [rp_hyp,a_hyp,e_hyp,deltadegree,deltav_pnorm,delta] = hyperbola(vinfminus,ip(k),mu_E);

rp_hyp_=[rp_hyp_, rp_hyp];
a_hyp_=[a_hyp_, a_hyp];
e_hyp_=[e_hyp_, e_hyp];
deltav_pnorm_=[deltav_pnorm_,deltav_pnorm];
delta_=[delta_,delta];
end

for k=1:length(ip)
[v_rot]= Rodrigues_vec(uset,vinfminus,delta_(k));
v_rot_ = [v_rot_,v_rot];
end

VPversor = cross(r_E,uset)./norm(cross(r_E,uset));
VP = sqrt(mu_S./norm(r_E)).*VPversor;

VPLUSnew=[];
for k=1:length(ip)
VPLUS= VP+v_rot_(:,k);
VMINUS = VP+vinfminus;
VPLUSnew = [VPLUSnew,VPLUS];
end
%% HELIO
figure()

y0_helio1=[];
y0_helio2= [];
for k=1:length(ip)
y01 = [r_E;VMINUS];
y02 = [r_E;VPLUSnew(:,k)];
y0_helio1= [y0_helio1,y01];
y0_helio2=[y0_helio2,y02];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
 
tspan1 = linspace(0,-365*24*3600/2,100);
tspan2=linspace(0,365*24*3600/2,100);

[TE, Y_helio1] = ode113( @(t,y) ode_2bp(t,y, mu_S), tspan1, y0_helio1(:,k), options );
[TM, Y_helio2] = ode113( @(t,y) ode_2bp(t,y, mu_S), tspan2, y0_helio2(:,k), options );
  
plot3( Y_helio1(:,1), Y_helio1(:,2), Y_helio1(:,3),'b' );
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
axis equal;
grid on;
hold on;
plot3( Y_helio2(:,1), Y_helio2(:,2), Y_helio2(:,3) );
hold on;
end
%% HYP
rpvector_=[];
for k=1:length(ip)
deltav_direction=v_rot_(:,k)-vinfminus;
rpvector=rp_hyp_(k).*(-deltav_direction./deltav_pnorm_(k));
rpvector_ = [rpvector_,rpvector];
end

vp_hyp_=[];
for k=1:length(ip)
vp_hyp = sqrt(norm(vinfminus).^2+2*mu_E./norm(rpvector_(:,k))).*(cross(rpvector_(:,k),uset)./norm(cross(rpvector_(:,k),uset)));
vp_hyp_ = [vp_hyp_,vp_hyp];
end

y0_hyp = [];
for k=1:length(ip)
y0 = [rpvector_(:,k);vp_hyp_(:,k)];
y0_hyp=[y0_hyp,y0];
end

figure()
for k=1:length(ip)

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
 
tspan1=linspace(0,-3000,100);
tspan2=linspace(0,3000,100);

[TE, Y_hyp1] = ode113( @(t,y) ode_2bp(t,y, mu_E), tspan1, y0_hyp(:,k), options );
[TE, Y_hyp2] = ode113( @(t,y) ode_2bp(t,y, mu_E), tspan2, y0_hyp(:,k), options );
 
plot3( Y_hyp1(:,1), Y_hyp1(:,2), Y_hyp1(:,3));
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
axis equal;
grid on;
hold on;
plot3( Y_hyp2(:,1), Y_hyp2(:,2), Y_hyp2(:,3));
plot3(rpvector_(1,k),rpvector_(2,k),rpvector_(3,k),'o');
 
Terra3d
hold on;
end

