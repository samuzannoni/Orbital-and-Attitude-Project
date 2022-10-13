clear all
clc
close all

Ixx=0.07;
Iyy=0.0504;
Izz=0.0109;

J=[Ixx 0 0; 0 Iyy 0; 0 0 Izz]; % inertia matrix
Jinv=inv(J); 

w0=[0.45 0.52 0.55]'; % initial conditions on angular velocity
A0=eye(3); % A0 needs to be otrthonormal because it represents the starting coordinates in an orthonormal reference frame

simulation=sim('lab4_dcm.slx');

t=simulation.tout;

w_hat=simulation.simout1;

A=simulation.simout2;

A_ortho=simulation.simout3;

A_check=simulation.simout4; % A_check is used to check the orthonormalisation: A_check=A_ortho*A_ortho'

bool=[];

for i=1:length(t)

    if (det(A_check(:,:,i))>(1-1e-4)) && (det(A_check(:,:,i))<(1+1e-4)) % 1e-4 is a toll used to check the orthonormality if A_check is not exactly equal to I
        
        bool(i)=0;
        
    else
        bool(i)=1;
    end
end

if norm(bool)==0
    disp('all matrixes A_ortho are orthonormal');
else
    disp('all qmatrixes A_ortho are not orthonormal');
end

