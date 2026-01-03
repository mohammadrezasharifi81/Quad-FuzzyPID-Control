% --- Script to Generate Path: Takeoff, Figure-Eight, Land ---

clear;
clc;

%% 1. Parameters
height = 5;         % meters
x_radius = 4;       % meters (width of the eight)
y_radius = 2;       % meters (height of the eight lobes)
sample_rate = 50;   % Hz
takeoff_time = 5;
figure_eight_time = 20;
land_time = 5;

%% 2. Time Vectors
t1 = (0:1/sample_rate:takeoff_time)';
t2 = (takeoff_time + 1/sample_rate : 1/sample_rate : takeoff_time + figure_eight_time)';
t3 = (t2(end) + 1/sample_rate : 1/sample_rate : t2(end) + land_time)';
time_vector = [t1; t2; t3];
angle = linspace(0, 2*pi, length(t2))';

%% 3. Coordinates
% Takeoff
x_takeoff = zeros(size(t1));
y_takeoff = zeros(size(t1));
z_takeoff = linspace(0, height, length(t1))';
% Figure-Eight
x_figure_eight = x_radius * sin(angle);
y_figure_eight = y_radius * sin(2*angle);
z_figure_eight = ones(size(t2)) * height;
% Land
x_land = ones(size(t3)) * x_figure_eight(end);
y_land = ones(size(t3)) * y_figure_eight(end);
z_land = linspace(height, 0, length(t3))';
% Combine
x_path = [x_takeoff; x_figure_eight; x_land];
y_path = [y_takeoff; y_figure_eight; y_land];
z_path = [z_takeoff; z_figure_eight; z_land];
psi_path = zeros(size(time_vector));

%% 4. Create Structure
path.x = timeseries(x_path, time_vector);
path.y = timeseries(y_path, time_vector);
path.z = timeseries(z_path, time_vector);
path.psi = timeseries(psi_path, time_vector);

%% 5. Save File
save('Path_Takeoff_FigureEight_Land.mat', 'path');
disp('File "Path_Takeoff_FigureEight_Land.mat" created successfully.');

%% 6. Plot
figure('Name', 'Generated Path: Figure-Eight');
plot3(x_path, y_path, z_path, 'b-', 'LineWidth', 2);
hold on;
plot3(x_path(1), y_path(1), z_path(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
plot3(x_path(end), y_path(end), z_path(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
hold off;
grid on; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Generated Path: Takeoff, Figure-Eight, Land'); axis equal; view(45, 25);