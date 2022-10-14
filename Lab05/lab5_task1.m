clear all
clc
close all

Ixx=0.07;
Iyy=0.0504;
Izz=0.0109;

J=[Ixx 0 0; 0 Iyy 0; 0 0 Izz]; % inertia matrix
Jinv=inv(J); 

w0=[0.45 0.52 0.55]'; % initial conditions on angular velocity

% starting set of Euler's angles: set them zero means body frame is coincident with the non inertial frame
phi0=0;
theta0=0;
psi0=0;

simulation=sim('lab5_EAs.slx');

t=simulation.tout;

wx=simulation.simout(:,1); % extrapolate the oemgas because i need to use them in phi, theta, psi
wy=simulation.simout(:,2);
wz=simulation.simout(:,3);

phi=simulation.EAs(:,1);
theta=simulation.EAs(:,2);
psi=simulation.EAs(:,3);

EA = [phi theta psi];

% A312=[cos(psi)cos(phi)-sin(psi)*sin(phi)*sin(theta) cos(psi)*sin(phi)+sin(psi)*cos(phi)*sin(theta) -sin(psi)*cos(theta);
%     -sin(phi)*cos(theta) cos(phi)*cos(theta) sin(theta);
%     sin(psi)*cos(phi)+cos(psi)*sin(phi)*sin(theta) sin(psi)*sin(phi)-cos(psi)*cos(phi)*sin(theta) cos(theta)*cos(psi)];

