close all
clc
clear

%% Task 1: general asymmetric spacecraft

Ixx=0.07;
Iyy=0.0504;
Izz=0.0109;

J=[Ixx 0 0; 0 Iyy 0; 0 0 Izz]; % inertia matrix
Jinv=inv(J); 

w0=[0.45 0.52 0.55]'; % initial conditions on angular velocity

simulation=sim('euler_equation.slx');
t=simulation.tout;
wx=simulation.simout(:, 1);
wy=simulation.simout(:, 2);
wz=simulation.simout(:, 3);

% plot the results (they must be equal to that of the scope)
figure(1)
plot(t, wx);
hold on;
grid on;
plot(t, wy);
plot(t, wz);
title('Asymmetrical Spacecraft - Simulation Results');
xlabel('time');
ylabel('angular velocity components');
legend('omega_x', 'omega_y', 'omega_z');

% %% Task 2: axial-symmetrical spacecraft
% 
% Ixx=0.0504; % Ixx = Iyy
% 
% % exact solution of the Euler equation:
% lambda=((Izz-Ixx)/Ixx)*w0(3);
% wx_ex=@(t) w0(1)*cos(lambda*t)-w0(2)*sin(lambda*t);
% wy_ex=@(t) w0(1)*sin(lambda*t)+w0(2)*cos(lambda*t);
% wz_ex=@(t) w0(3)+t*0;
% 
% % plot the exact solution
% figure(2);
% plot(t, wx_ex(t));
% hold on;
% grid on;
% plot(t, wy_ex(t));
% plot(t, wz_ex(t));
% title('Axial-symmetrical Spacecraft - Analitical Solutions');
% xlabel('time');
% ylabel('angular velocity components');
% legend('omega_x', 'omega_y', 'omega_z');
% 
% % plot the Absolute Difference between Simulation Results and Exact Solutions
% figure(3);
% plot(t, abs(wx-wx_ex(t)));
% hold on;
% grid on;
% plot(t, abs(wy-wy_ex(t)));
% plot(t, abs(wz-wz_ex(t)));
% title('Axial-symmetrical Spacecraft - Absolute Error');
% xlabel('time');
% ylabel('abs(w_sim - w_ex)');
% legend('delta_omega_x', 'delta_omega_y', 'delta_omega_z');
% 
% % max error made with the simulation
% err_max_wx=max(abs(wx-wx_ex(t)));
% err_max_wy=max(abs(wy-wy_ex(t)));
% err_max_wz=max(abs(wz-wz_ex(t)));
% err_max=[err_max_wx err_max_wy err_max_wz];
% 
% % symmetric exact solution does provide a good approximation of the
% % asymmetric case????

%% Task 3: stability

% Ixx=0.01; % X-axis is minor axis
% Iyy=0.05; % Y-axis is intermediate inertial axis
% Izz=0.07; % Z-axis is major inertial axis
% 
% J=[Ixx 0 0; 0 Iyy 0; 0 0 Izz]; % inertia matrix
% Jinv=inv(J);
% 
% w0_matrix=2*pi*eye(3);
% w_dot_0_matrix=zeros(3);
% 
%  simulation=sim('euler_equation.slx');
%  tspan=simulation.tout';
% wx=simulation.simout(:, 1);
% wy=simulation.simout(:, 2);
% wz=simulation.simout(:, 3);
% 
% % verify that angular momentum and kinetic energy are constant
% hh=[Ixx*wx Iyy*wy Izz*wz]; 
% 
% h=zeros(size(hh,1),1);
% for i=1:size(hh,1)
%     h(i)=norm(hh(i,:)); % 
% end
% 
% T=0.5*(((hh(:,1).^2)./Ixx+(hh(:,2)).^2)./Iyy+(hh(:,3).^2)./Izz);

% kinetic_energy_ellipsoid=@(wx,wy,wz) (wx.^2)./((2*T)./Ixx)+(wy.^2)./((2*T)./Iyy)+(wz.^2)./((2*T)./Izz)-1;

% angular_momentum_ellipsoid=@(wx,wy,wz) (wx.^2)./((h.^2)./Ixx)+(wy.^2)./((h.^2)./Iyy)+(wz.^2)./((h.^2)./Izz)-1;

% figure(4);
% [Xk, Yk, Zk]=ellipsoid(w0(1),w0(2),w0(3),2*T(1)/Ixx,2*T(1)/Iyy,2*T(1)/Izz); % kinetic_energy_ellipsoid
% surf(Xk,Yk,Zk);
% axis equal;
% hold on;
% [Xh, Yh, Zh]=ellipsoid(w0(1),w0(2),w0(3),(h(1)^2)./(Ixx^2),(h(1)^2)./(Iyy^2),(h(1)^2)./(Izz^2)); % angular_momentum_ellipsoid
% surf(Xh,Yh,Zh);
% title('Ellipsoid Intersection');


% switch n='1, 2, 3'
% 
%     case 1
%     n=2;
%         w0=w0_matrix(n,:)';
%         K=((Ixx-Izz)*(Ixx-Iyy)*w0(n)^2)/(Iyy*Izz);
%         w_dot_0=w_dot_0_matrix(n,:)';
%         y0=[w0; w_dot_0];
%         [ T, Y] = ode113( @(t,y) stability(t,y,K), tspan, y0, options);
%         figure(4)
%         plot(T, Y(:,1));
%         
% %     case 2
% %         w0=w0_matrix(w0(n,:));
