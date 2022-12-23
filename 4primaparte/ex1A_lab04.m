%Design a flyby around the Earth for fixed impact parameter and different locations of 
%the incoming asymptote
%% HYP planetocentric
%1.A
mu_E = astroConstants(13);
mu_S = astroConstants(4);
AU = astroConstants(2);
vinfminus=[15.1; 0; 0];
ip = 9200; % impact parameter
r_E=[1; 0; 0]*AU;

% Charachterization of the hyperbola
[rp_hyp,a_hyp,e_hyp,deltadegree,deltav_pnorm,delta] = hyperbola(vinfminus,ip,mu_E);

% v_inf_plus 
ubehind= [0; 0; 1];
ufront= [0; 0; -1]; % la terra gira counter clockwise
ubelow= [0; -1; 0];
u = [ubehind,ufront,ubelow];

[v_rotated]= Rodrigues(u,vinfminus,delta);

vinfplus_behind=v_rotated(:,1); 
vinfplus_front=v_rotated(:,2);
vinfplus_below=v_rotated(:,3);

VPversor = cross(r_E,ufront)./norm(cross(r_E,ufront));
VP = sqrt(mu_S./norm(r_E)).*VPversor;

% VPLUS
VPLUS = [];
for i= 1:length(u)
VPLUS_ = VP+v_rotated(:,i);
VPLUS= [VPLUS,VPLUS_];
end

VPLUS_behind=VPLUS(:,1); 
VPLUS_front=VPLUS(:,2);
VPLUS_below=VPLUS(:,3);

VMINUS = VP+vinfminus;

%Plot the trajectories
y0_hyp=[];
rpvector_=[];

for i=1:length(u)
deltav_direction=v_rotated(:,i)-vinfminus;
rpvector=rp_hyp.*(-deltav_direction./deltav_pnorm);
vp_hyp = sqrt(norm(vinfminus).^2+2*mu_E./norm(rpvector)).*(cross(rpvector,u(:,i))./norm(cross(rpvector,u(:,i))));
y0 = [rpvector;vp_hyp];
y0_hyp=[y0_hyp,y0];
rpvector_ = [rpvector_,rpvector];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

tspan1=linspace(0,3000,100);
tspan2=linspace(0,-3000,100);

[TE, Y_hyp1] = ode113( @(t,y) ode_2bp(t,y, mu_E), tspan1, y0_hyp(:,i), options );
[TE, Y_hyp2] = ode113( @(t,y) ode_2bp(t,y, mu_E), tspan2, y0_hyp(:,i), options );

figure()
plot3( Y_hyp1(:,1), Y_hyp1(:,2), Y_hyp1(:,3),'b','LineWidth',2 );
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
axis equal;
grid on;
hold on;
plot3( Y_hyp2(:,1), Y_hyp2(:,2), Y_hyp2(:,3),'r','LineWidth',2 );
plot3(rpvector_(1,i),rpvector_(2,i),rpvector_(3,i),'og','Linewidth',2,'MarkerSize',8);

Terra3d

end
%% HELIO

y0_helio1=[];
y0_helio2 = [];

for i=1:length(u)
y01 = [r_E;VMINUS];
y02 = [r_E;VPLUS(:,i)];
y0_helio1=[y0_helio1,y01];
y0_helio2= [y0_helio2,y02];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

tspan1 = linspace(0,-365*24*3600/2,100);
tspan2=linspace(0,365*24*3600/2,100);

% Perform the integration
[TE, Y_helio1] = ode113( @(t,y) ode_2bp(t,y, mu_S), tspan1, y0_helio1(:,i), options );
[TM, Y_helio2] = ode113( @(t,y) ode_2bp(t,y, mu_S), tspan2, y0_helio2(:,i), options );


%Plot
figure
plot3( Y_helio1(:,1), Y_helio1(:,2), Y_helio1(:,3),'b' );
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
axis equal;
grid on;
hold on;
plot3( Y_helio2(:,1), Y_helio2(:,2), Y_helio2(:,3),'r' );
end
