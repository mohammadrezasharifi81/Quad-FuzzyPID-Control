function [] = quadPlots(yout,tout)
% Clear all previous figures and the command window
close all; clc;

% Assign input data to local variables
A = yout;
T = tout;
t = T;

% Extract state variables from the input matrix 'A'
P = A(:,1); 
Q = A(:,2);
R = A(:,3);
Phi   = A(:,4); 
Theta = A(:,5);
Psi   = A(:,6);
U = A(:,7);
V = A(:,8);
W = A(:,9);
X = A(:,10);
Y = A(:,11);
Z = A(:,12);

% Extract motor speeds (RPM) and motor commands (Throttle %)
w1 = A(:,13);
w2 = A(:,14);
w3 = A(:,15);
w4 = A(:,16);
mc1 = A(:,17);
mc2 = A(:,18);
mc3 = A(:,19);
mc4 = A(:,20);

% Extract commanded (reference) values
Phi_cmd = A(:,21);
Theta_cmd = A(:,22);
Psi_cmd = A(:,23);
Z_cmd   = A(:,24);

% Check if the simulation includes position commands for X and Y
[~,column] = size(yout);
if column==26
    PC = true;
    X_cmd = A(:,25);
    Y_cmd = A(:,26);
else
    PC = false;
end

% --- Figure 1: Position and Attitude Plots ---
figure('Name', 'Quadcopter State Variables');

% Plot for X position
subplot(2,3,1)
plot(T,X,'b')
if (PC==true)
    hold on
    plot(T,X_cmd,'k--')
    hold off
end
xlabel('Time (s)')
ylabel('X Position (m)')
xlim([min(t) max(t)])
title('X Position')
grid on

% Plot for Y position
subplot(2,3,2)
plot(T,Y,'r')
if (PC==true)
    hold on
    plot(T,Y_cmd,'k--')
    hold off
end
xlabel('Time (s)')
ylabel('Y Position (m)')
xlim([min(t) max(t)])
title('Y Position')
grid on

% Plot for Z position (Altitude)
subplot(2,3,3)
plot(T,Z,'g')
hold on
plot(T,Z_cmd,'k--')
hold off
xlabel('Time (s)')
ylabel('Z Position (m)')
xlim([min(t) max(t)])
title('Z Position (Altitude)')
grid on

% Plot for Phi angle (Roll)
subplot(2,3,4)
plot(T,Phi*180/pi,'b')
hold on
plot(T,Phi_cmd*180/pi,'k--')
hold off
xlabel('Time (s)')
ylabel('Angle (deg)')
xlim([min(t) max(t)])
title('Phi (Roll)')
grid on
 
% Plot for Theta angle (Pitch)
subplot(2,3,5)
plot(T,Theta*180/pi,'r')
hold on
plot(T,Theta_cmd*180/pi,'k--')
hold off
xlabel('Time (s)')
ylabel('Angle (deg)')
xlim([min(t) max(t)])
title('Theta (Pitch)')
grid on

% Plot for Psi angle (Yaw)
subplot(2,3,6)
plot(T,Psi*180/pi,'g')
hold on
plot(T,Psi_cmd*180/pi,'k--')
hold off
xlabel('Time (s)')
ylabel('Angle (deg)')
xlim([min(t) max(t)])
title('Psi (Yaw)')
grid on

% --- Figure 2: Motor RPM and Control Effort ---
figure('Name', 'Motor Performance');

% Motor 1
subplot(4,1,1)
[AX,H1,H2] = plotyy(T,mc1,T,w1);
xlabel('Time (s)')
ylabel(AX(1),'Motor Throttle Command (%)')
ylabel(AX(2),'Motor Speed (RPM)')
xlim(AX(1),[min(T) max(T)])
xlim(AX(2),[min(T) max(T)])
title('Motor 1')
grid on

% Motor 2
subplot(4,1,2)
[AX,H1,H2] = plotyy(T,mc2,T,w2);
xlabel('Time (s)')
ylabel(AX(1),'Motor Throttle Command (%)')
ylabel(AX(2),'Motor Speed (RPM)')
xlim(AX(1),[min(T) max(T)])
xlim(AX(2),[min(T) max(T)])
title('Motor 2')
grid on

% Motor 3
subplot(4,1,3)
[AX,H1,H2] = plotyy(T,mc3,T,w3);
xlabel('Time (s)')
ylabel(AX(1),'Motor Throttle Command (%)')
ylabel(AX(2),'Motor Speed (RPM)')
xlim(AX(1),[min(T) max(T)])
xlim(AX(2),[min(T) max(T)])
title('Motor 3')
grid on

% Motor 4
subplot(4,1,4)
[AX,H1,H2] = plotyy(T,mc4,T,w4);
xlabel('Time (s)')
ylabel(AX(1),'Motor Throttle Command (%)')
ylabel(AX(2),'Motor Speed (RPM)')
xlim(AX(1),[min(T) max(T)])
xlim(AX(2),[min(T) max(T)])
title('Motor 4')
grid on

end