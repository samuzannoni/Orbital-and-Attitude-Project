close all
clc
clear

%% Task 2: axial-symmetrical spacecraft

Ixx=0.0504; % Ixx = Iyy
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

% exact solution of the Euler equation:
lambda=((Izz-Ixx)/Ixx)*w0(3);
wx_ex=@(t) w0(1)*cos(lambda*t)-w0(2)*sin(lambda*t);
wy_ex=@(t) w0(1)*sin(lambda*t)+w0(2)*cos(lambda*t);
wz_ex=@(t) w0(3)+t*0;

% plot the exact solution
figure(2);
plot(t, wx_ex(t));
hold on;
grid on;
plot(t, wy_ex(t));
plot(t, wz_ex(t));
title('Axial-symmetrical Spacecraft - Analitical Solutions');
xlabel('time');
ylabel('angular velocity components');
legend('omega_x', 'omega_y', 'omega_z');

% plot the Absolute Difference between Simulation Results and Exact Solutions
figure(3);
plot(t, abs(wx-wx_ex(t)));
hold on;
grid on;
plot(t, abs(wy-wy_ex(t)));
plot(t, abs(wz-wz_ex(t)));
title('Axial-symmetrical Spacecraft - Absolute Error');
xlabel('time');
ylabel('abs(w_sim - w_ex)');
legend('delta_omega_x', 'delta_omega_y', 'delta_omega_z');

% max error made with the simulation
err_max_wx=max(abs(wx-wx_ex(t)));
err_max_wy=max(abs(wy-wy_ex(t)));
err_max_wz=max(abs(wz-wz_ex(t)));
err_max=[err_max_wx err_max_wy err_max_wz];

% symmetric exact solution does provide a good approximation of the
% asymmetric case????