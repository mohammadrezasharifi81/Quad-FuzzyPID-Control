clear; clc; close all;

%% ========================================================================
%  1. SETUP AND DATA LOADING
%  ========================================================================
% Load the simulation results file
fileName = 'simResults.mat';
if ~exist(fileName, 'file')
    error('File "%s" not found. Make sure it is in the current folder.', fileName);
end

load(fileName);
fprintf('Loaded simulation results from: %s\n', fileName);

% Define the time vector for plotting
time = tout;

%% ========================================================================
%  2. DATA EXTRACTION AND PREPARATION
%  ========================================================================
% Determine the simulation type based on the number of output columns
[~, numColumns] = size(yout);
isPositionControl = (numColumns == 26);

% Extract Actual (Simulated) States
phi_actual   = yout(:, 4);   % Roll [rad]
theta_actual = yout(:, 5);   % Pitch [rad]
psi_actual   = yout(:, 6);   % Yaw [rad]
x_actual     = yout(:, 10);  % X [m]
y_actual     = yout(:, 11);  % Y [m]
z_actual     = yout(:, 12);  % Z [m]

% Extract Desired (Reference/Command) States
phi_ref_cmd   = yout(:, 21); % Commanded Roll [rad]
theta_ref_cmd = yout(:, 22); % Commanded Pitch [rad]
psi_ref_cmd   = yout(:, 23); % Commanded Yaw [rad]

% Convert attitude angles from radians to degrees (for readability)
phi_actual   = rad2deg(phi_actual);
theta_actual = rad2deg(theta_actual);
psi_actual   = rad2deg(psi_actual);

phi_ref_cmd   = rad2deg(phi_ref_cmd);
theta_ref_cmd = rad2deg(theta_ref_cmd);
psi_ref_cmd   = rad2deg(psi_ref_cmd);

% For position control, interpolate the high-level path to the simulation time
if isPositionControl
    x_desired = interp1(path.x.Time, squeeze(path.x.Data), time);
    y_desired = interp1(path.y.Time, squeeze(path.y.Data), time);
    z_data_raw = squeeze(path.z.Data);
    if numel(z_data_raw) == 1
        z_desired = repmat(z_data_raw, size(x_desired));
    else
        z_desired = interp1(path.z.Time, z_data_raw, time);
    end
end

%% ========================================================================
%  3. 3D FLIGHT ANIMATION  -->  SAVE LOWER-QUALITY VIDEO
%  ========================================================================
disp('Starting animation (video capture)...');

% ---------- Appearance / geometry ----------
arm_length = 0.25;
body_size  = 0.05;
quad_arms  = [-arm_length, 0, 0;
               arm_length, 0, 0;
               0, -arm_length, 0;
               0,  arm_length, 0];
quad_body  = [ body_size,  body_size, 0;
              -body_size,  body_size, 0;
              -body_size, -body_size, 0;
               body_size, -body_size, 0;
               body_size,  body_size, 0];

% ---------- Video settings (reduced resolution & quality) ----------
outName = 'quad_anim_540p.mp4'; % Output file name
FPS     = 10;                   % Video frame rate
QUALITY = 25;                   % 1..100 (lower = smaller file / lower quality)
OUT_RES = [960 540];            % Output resolution [width height] in pixels

% ---------- Figure (smaller window, matches output resolution) ----------
fig_anim = figure('Name','Quadcopter Animation', ...
                  'Color','w', ...
                  'Position',[100 100 OUT_RES]); % Window size ~ video resolution
set(fig_anim,'Renderer','opengl'); % Better 3D rendering

% ---------- Axes prepared once (faster than creating subplots each frame) ----------
tiledlayout(fig_anim,1,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal'); view(ax1,30,25);
xlabel(ax1,'X (m)'); ylabel(ax1,'Y (m)'); zlabel(ax1,'Z (m)');
xlim(ax1,[-0.5 0.5]); ylim(ax1,[-0.5 0.5]); zlim(ax1,[-0.5 0.5]);
title(ax1,'Attitude');

ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); axis(ax2,'equal'); view(ax2,30,25);
xlabel(ax2,'X Position (m)'); ylabel(ax2,'Y Position (m)'); zlabel(ax2,'Z Position (m)');
title(ax2,'Position');
if isPositionControl
    plot3(ax2, x_desired, y_desired, z_desired, 'r--', 'LineWidth', 1.2); % Reference path
end
hPath   = plot3(ax2, NaN,NaN,NaN, 'b-',  'LineWidth', 2.0);                % Traveled path
hMarker = plot3(ax2, NaN,NaN,NaN, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);

% Body/arm graphics objects to be updated each frame
hArm12 = plot3(ax1, NaN,NaN,NaN, 'k-', 'LineWidth', 2.2);
hArm34 = plot3(ax1, NaN,NaN,NaN, 'k-', 'LineWidth', 2.2);
hAxis  = plot3(ax1, NaN,NaN,NaN, 'r-', 'LineWidth', 3.0);
hBody  = fill3(ax1, NaN,NaN,NaN, [0.4 0.4 0.4]);

% ---------- VideoWriter ----------
try
    v = VideoWriter(outName,'MPEG-4');  % H.264 (mp4)
catch
    warning('MPEG-4 profile not available. Falling back to Motion JPEG AVI.');
    outName = 'quad_anim_540p.avi';
    v = VideoWriter(outName,'Motion JPEG AVI');
end
v.FrameRate = FPS;
if strcmp(v.FileFormat,'MPEG-4'), v.Quality = QUALITY; end
open(v);

% ---------- Target frame timing (keep near-real-time speed) ----------
dt_target = 1/FPS;
if numel(time) > 1
    t0 = time(1); tEnd = time(end);
else
    t0 = 0;       tEnd = length(phi_actual)/50; % Fallback guess
end
t_frames = t0:dt_target:tEnd;

% Index in the simulation data for each target frame
k_idx = arrayfun(@(tf) find(time>=tf,1,'first'), t_frames);
k_idx(isnan(k_idx)) = numel(time);

% ---------- Animation loop + capture ----------
for ii = 1:numel(k_idx)
    k = k_idx(ii);

    % Angles (deg->rad) and rotation matrix (Z-Y-X convention)
    phi   = deg2rad(phi_actual(k));
    theta = deg2rad(theta_actual(k));
    psi   = deg2rad(psi_actual(k));
    Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    R  = Rz*Ry*Rx;

    % Rotate body/arms
    rotated_arms = (R * quad_arms')'; % 4x3
    rotated_body = (R * quad_body')'; % 5x3

    % Update graphics efficiently
    set(hArm12,'XData',rotated_arms([1,2],1),'YData',rotated_arms([1,2],2),'ZData',rotated_arms([1,2],3));
    set(hArm34,'XData',rotated_arms([3,4],1),'YData',rotated_arms([3,4],2),'ZData',rotated_arms([3,4],3));
    set(hAxis, 'XData',[0 rotated_arms(2,1)], 'YData',[0 rotated_arms(2,2)], 'ZData',[0 rotated_arms(2,3)]);
    set(hBody, 'XData',rotated_body(:,1),'YData',rotated_body(:,2),'ZData',rotated_body(:,3));

    title(ax1, sprintf('Attitude (Roll: %.1f°, Pitch: %.1f°, Yaw: %.1f°)', ...
        phi_actual(k), theta_actual(k), psi_actual(k)));

    % Path up to index k
    set(hPath,  'XData',x_actual(1:k),'YData',y_actual(1:k),'ZData',z_actual(1:k));
    set(hMarker,'XData',x_actual(k),'YData',y_actual(k),'ZData',z_actual(k));
    title(ax2, sprintf('Position (Time: %.2f s)', time(k)));

    drawnow limitrate nocallbacks;  % Fast rendering
    F = getframe(fig_anim);         % Capture frame (at OUT_RES)
    writeVideo(v, F);
end
close(v);
fprintf('Animation saved to: %s (FPS=%d, Quality=%d)\n', outName, FPS, QUALITY);

%% ========================================================================
%  4. STATIC PLOTS FOR PERFORMANCE ANALYSIS (smaller window)
%  ========================================================================
disp('Plotting 6-DOF tracking performance...');
figure('Name', 'Position and Angles Tracking', 'Position', [160, 160, 700, 500]);

% --- Position Plots ---
subplot(3, 2, 1);
if isPositionControl, plot(time, x_desired, 'b-', 'LineWidth', 1.2); hold on; end
plot(time, x_actual, 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('X (m)');
title('Position Tracking'); if isPositionControl, legend('Reference', 'Actual'); end

subplot(3, 2, 3);
if isPositionControl, plot(time, y_desired, 'b-', 'LineWidth', 1.2); hold on; end
plot(time, y_actual, 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('Y (m)');
if isPositionControl, legend('Reference', 'Actual'); end

subplot(3, 2, 5);
if isPositionControl, plot(time, z_desired, 'b-', 'LineWidth', 1.2); hold on; end
plot(time, z_actual, 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('Z (m)');
if isPositionControl, legend('Reference', 'Actual'); end

% --- Angle Plots (degrees) ---
subplot(3, 2, 2);
plot(time, phi_ref_cmd, 'b-', 'LineWidth', 1.2); hold on;
plot(time, phi_actual, 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('\phi (Roll) [deg]');
title('Angles Tracking'); legend('Command', 'Actual');

subplot(3, 2, 4);
plot(time, theta_ref_cmd, 'b-', 'LineWidth', 1.2); hold on;
plot(time, theta_actual, 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('\theta (Pitch) [deg]');
legend('Command', 'Actual');

subplot(3, 2, 6);
plot(time, psi_ref_cmd, 'b-', 'LineWidth', 1.2); hold on;
plot(time, psi_actual, 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('\psi (Yaw) [deg]');
legend('Command', 'Actual');

sgtitle('Position and Angles Tracking Performance');
