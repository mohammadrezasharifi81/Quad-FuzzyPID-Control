% --- Script to Generate Path: Diamond with Return and Landing at Origin ---

clear;
clc;

%% 1. Parameters
height = 3;         % meters
points_per_leg = 100; % Number of points for each segment
takeoff_time = 4;
return_time = 3;    % --- NEW: Time to return to center ---
land_time = 4;      % Time for landing phase
sample_rate = 50;

% Diamond vertices [x, y, z]
v1 = [0, 1, height];
v2 = [-1, 0, height];
v3 = [0, -1, height];
v4 = [1, 0, height];
vertices = [v1; v2; v3; v4; v1]; % Path ends back at v1

%% 2. Coordinates & Time
% --- Phase 1: Takeoff ---
t_takeoff = (0:1/sample_rate:takeoff_time)';
x_takeoff = zeros(size(t_takeoff));
y_takeoff = zeros(size(t_takeoff));
z_takeoff = linspace(0, height, length(t_takeoff))';

% --- Phase 2: Diamond Legs ---
x_diamond = []; y_diamond = []; z_diamond = [];
for i = 1:size(vertices, 1)-1
    x_leg = linspace(vertices(i,1), vertices(i+1,1), points_per_leg)';
    y_leg = linspace(vertices(i,2), vertices(i+1,2), points_per_leg)';
    z_leg = ones(points_per_leg, 1) * height;
    
    x_diamond = [x_diamond; x_leg(2:end)];
    y_diamond = [y_diamond; y_leg(2:end)];
    z_diamond = [z_diamond; z_leg(2:end)];
end
total_diamond_points = length(x_diamond);
diamond_duration = (total_diamond_points / sample_rate);
t_diamond = linspace(t_takeoff(end) + 1/sample_rate, t_takeoff(end) + diamond_duration, total_diamond_points)';

% --- Phase 3: Return to Center (NEW) ---
t_return = (t_diamond(end) + 1/sample_rate : 1/sample_rate : t_diamond(end) + return_time)';
x_return = linspace(x_diamond(end), 0, length(t_return))'; % Fly from diamond's end X to 0
y_return = linspace(y_diamond(end), 0, length(t_return))'; % Fly from diamond's end Y to 0
z_return = ones(size(t_return)) * height;                 % Maintain altitude

% --- Phase 4: Land at Origin (Modified) ---
t_land = (t_return(end) + 1/sample_rate : 1/sample_rate : t_return(end) + land_time)';
x_land = zeros(size(t_land)); % Land at X=0
y_land = zeros(size(t_land)); % Land at Y=0
z_land = linspace(height, 0, length(t_land))'; % Descend to Z=0

% --- Combine all phases ---
x_path = [x_takeoff; x_diamond; x_return; x_land];
y_path = [y_takeoff; y_diamond; y_return; y_land];
z_path = [z_takeoff; z_diamond; z_return; z_land];
time_vector = [t_takeoff; t_diamond; t_return; t_land];
psi_path = zeros(size(time_vector));

%% 3. Create Structure
path.x = timeseries(x_path, time_vector);
path.y = timeseries(y_path, time_vector);
path.z = timeseries(z_path, time_vector);
path.psi = timeseries(psi_path, time_vector);

%% 4. Save File
save('Path_Diamond_Return_Land.mat', 'path');
disp('File "Path_Diamond_Return_Land.mat" created successfully.');

%% 5. Plot
figure('Name', 'Generated Path: Diamond with Return and Landing');
plot3(x_path, y_path, z_path, 'b-', 'LineWidth', 2);
hold on;
plot3(x_path(1), y_path(1), z_path(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
plot3(x_path(end), y_path(end), z_path(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
hold off;
grid on; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Generated Path: Diamond with Return and Landing'); axis equal; view(30, 20);