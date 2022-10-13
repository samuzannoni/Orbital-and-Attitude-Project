clear all
clc
close all

Ixx=0.07;
Iyy=0.0504;
Izz=0.0109;

J=[Ixx 0 0; 0 Iyy 0; 0 0 Izz]; % inertia matrix
Jinv=inv(J); 

w0=[0.45 0.52 0.55]'; % initial conditions on angular velocity
e=[1; 0; 0];
phi=30; %deg
phi=deg2rad(phi);
q0=[e(1)*sin(phi/2); e(2)*sin(phi/2); e(3)*sin(phi/2); cos(phi/2)];

simulation=sim('lab4_quaternions.slx');

t=simulation.tout;

omega=simulation.OMEGA;

q=simulation.q;

mag_q=simulation.check;

bool=[];

for i=1:length(t)

    if (mag_q(i)>(1-1e-3))&&(mag_q(i)<(1+1e-3))
        
        bool(i)=0;
        
    else

        bool(i)=1;

    end
end

if norm(bool)==0
    disp('all quaternions are normal');
else
    disp('all quaternions are not normal');
end

