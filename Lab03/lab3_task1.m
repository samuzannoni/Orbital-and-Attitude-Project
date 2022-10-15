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
