close all
clc
clear all

%% Task 3: stability

Ixx=0.01; % X-axis is minor inertial axis
Iyy=0.05; % Y-axis is intermediate inertial axis
Izz=0.07; % Z-axis is major inertial axis

J=[Ixx 0 0; 0 Iyy 0; 0 0 Izz]; % inertia matrix
Jinv=inv(J);

w0_matrix=2*pi*eye(3);

    %ROTATION AROUND MINOR AXIS

        w0=w0_matrix(1,:)';
        simulation=sim('euler_equation.slx');
        t=simulation.tout;
        wx=simulation.simout(:, 1);
        wy=simulation.simout(:, 2);
        wz=simulation.simout(:, 3);
    
        figure(1);
        plot(t, wx);
        hold on;
        grid on;
        plot(t, wy);
        plot(t, wz);
        legend('omega_x','omega_y','omega_z');
        title('ROTATION AROUND MINOR AXIS');
        xlabel('time [s]'); ylabel('omega [rad/s]');

     %ROTATION AROUND INTERMEDIATE AXIS

        w0=w0_matrix(2,:)';
        simulation=sim('euler_equation.slx');
        t=simulation.tout;
        wx=simulation.simout(:, 1);
        wy=simulation.simout(:, 2);
        wz=simulation.simout(:, 3);
    
        figure(2);
        plot(t, wx);
        hold on;
        grid on;
        plot(t, wy);
        plot(t, wz);
        legend('omega_x','omega_y','omega_z');
        title('ROTATION AROUND INTERMEDIATE AXIS');
        xlabel('time [s]'); ylabel('omega [rad/s]');

    %ROTATION AROUND MAJOR AXIS

        w0=w0_matrix(3,:)';
        simulation=sim('euler_equation.slx');
        t=simulation.tout;
        wx=simulation.simout(:, 1);
        wy=simulation.simout(:, 2);
        wz=simulation.simout(:, 3);
    
        figure(3);
        plot(t, wx);
        hold on;
        grid on;
        plot(t, wy);
        plot(t, wz);
        legend('omega_x','omega_y','omega_z');
        title('ROTATION AROUND MAJOR AXIS');
        xlabel('time [s]'); ylabel('omega [rad/s]');
