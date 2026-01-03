% --- Script to Generate a Path: Takeoff, Circle, Return, and Land at Origin ---
%% 1. Define Path Parameters
height = 2;         % meters
radius = 1.5;       % meters
sample_rate = 50;   % points per second (Hz)

takeoff_time = 4;   % seconds to reach altitude
circle_time = 10;   % seconds to complete the circle
return_time = 3;    % --- NEW: seconds to return to the center ---
land_time = 4;      % seconds to land

%% 2. Generate Time and Angle Vectors
% Create time vectors for each of the four phases
t1 = (0:1/sample_rate:takeoff_time)';
t2 = (t1(end) + 1/sample_rate : 1/sample_rate : t1(end) + circle_time)';
t3 = (t2(end) + 1/sample_rate : 1/sample_rate : t2(end) + return_time)'; % --- NEW: Return time ---
t4 = (t3(end) + 1/sample_rate : 1/sample_rate : t3(end) + land_time)'; % --- NEW: Landing time ---

% Full time vector
time_vector = [t1; t2; t3; t4];

% Angle vector for the circle
angle = linspace(0, 2*pi, length(t2))';

%% 3. Generate Coordinates for Each Phase
% --- Phase 1: Takeoff ---
x_takeoff = zeros(size(t1));
y_takeoff = zeros(size(t1));
z_takeoff = linspace(0, height, length(t1))';

% --- Phase 2: Circle ---
x_circle = radius * cos(angle);
y_circle = radius * sin(angle);
z_circle = ones(size(t2)) * height;

% --- Phase 3: Return to Center (NEW) ---
x_return = linspace(x_circle(end), 0, length(t3))'; % Fly from circle's end X to 0
y_return = linspace(y_circle(end), 0, length(t3))'; % Fly from circle's end Y to 0
z_return = ones(size(t3)) * height;                 % Maintain altitude

% --- Phase 4: Land at Origin (Modified) ---
x_land = zeros(size(t4)); % Land at X=0
y_land = zeros(size(t4)); % Land at Y=0
z_land = linspace(height, 0, length(t4))'; % Descend to Z=0

% --- Combine All Four Phases ---
x_path = [x_takeoff; x_circle; x_return; x_land];
y_path = [y_takeoff; y_circle; y_return; y_land];
z_path = [z_takeoff; z_circle; z_return; z_land];
psi_path = zeros(size(time_vector));

%% 4. Create the Path Structure with Timeseries Objects
path.x = timeseries(x_path, time_vector);
path.y = timeseries(y_path, time_vector);
path.z = timeseries(z_path, time_vector);
path.psi = timeseries(psi_path, time_vector);

%% 5. Save the Path to a new .mat File
save('Path_Takeoff_Circle.mat', 'path');
disp('File "Path_Takeoff_Circle.mat" created successfully.');

%% 6. Plot the Generated Path
try
    figure('Name', 'Generated 3D Path: Takeoff, Circle, Return, Land');
    plot3(x_path, y_path, z_path, 'b-', 'LineWidth', 2);
    
    hold on;
    plot3(x_path(1), y_path(1), z_path(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10); % Start
    plot3(x_path(end), y_path(end), z_path(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10); % End
    hold off;

    grid on;
    xlabel('X-axis (m)');
    ylabel('Y-axis (m)');
    zlabel('Z-axis (m)');
    title('Flight Trajectory: Takeoff, Circle, Return, Land');
    legend('Path', 'Start', 'End', 'Location', 'best');
    axis equal;
    view(30, 20);

    disp('3D plot generated successfully.');
catch ME
    warning('Could not generate 3D plot. Error: %s', ME.message);
end