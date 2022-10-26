clear all
clc
close all

Ixx=0.07;
Iyy=0.0504;
Izz=0.0109;

J=[Ixx 0 0; 0 Iyy 0; 0 0 Izz]; % inertia matrix
Jinv=inv(J);

C=0.2; % 2*pi<=C<=0.2
w0=[C 0.1 0.1]'; % initial conditions on angular velocity [rad/s]

simulation=sim('lab4_dcm.slx');

t=simulation.tout;

wx=simulation.omegas(:,1); % extrapolate the oemgas because i need to use them in phi, theta, psi
wy=simulation.omegas(:,2);
wz=simulation.omegas(:,3);

h_body=simulation.h_body;

v_test=[0 1 0]';
