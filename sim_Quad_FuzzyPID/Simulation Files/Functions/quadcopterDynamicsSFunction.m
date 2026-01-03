% ==========================
% Quadcopter Dynamics S-Function
% Author: [Mohammadreza Sharifi]
% Description:
%   Nonlinear dynamic model of a quadcopter implemented as a Simulink
%   Level-2 MATLAB S-Function. Includes 12 continuous states:
%   [P, Q, R, Phi, Theta, Psi, U, V, W, X, Y, Z]
% ==========================

function quadcopterDynamicsSFunction(block)
    setup(block);
end

% ==========================
% Initialization and Setup
% ==========================
function setup(block)

    % Define number of input and output ports
    block.NumInputPorts  = 5;
    block.NumOutputPorts = 12;

    % Configure input ports (4 motor speeds + 1 disturbance vector)
    for i = 1:4 
        block.InputPort(i).Dimensions        = 1;
        block.InputPort(i).DirectFeedthrough = false;
        block.InputPort(i).SamplingMode      = 'Sample';
    end

    block.InputPort(5).Dimensions        = 6; % [Dist_tau; Dist_F]
    block.InputPort(5).DirectFeedthrough = false;
    block.InputPort(5).SamplingMode      = 'Sample';

    % Configure output ports (12 states)
    for i = 1:12
        block.OutputPort(i).Dimensions   = 1;
        block.OutputPort(i).SamplingMode = 'Sample';
    end

    % Number of dialog parameters (quad structure and initial conditions)
    block.NumDialogPrms     = 2;
    block.NumContStates     = 12;

    block.SetAccelRunOnTLC(false);
    block.SimStateCompliance = 'DefaultSimState';

    % Register block methods
    block.RegBlockMethod('CheckParameters', @CheckPrms);
    block.RegBlockMethod('InitializeConditions', @InitializeConditions);
    block.RegBlockMethod('Outputs', @Outputs);
    block.RegBlockMethod('Derivatives', @Derivatives);
end

% ==========================
% Parameter Validation
% ==========================
function CheckPrms(block)
    quad = block.DialogPrm(1).Data;
    IC   = block.DialogPrm(2).Data;
end

% ==========================
% Initial Conditions Setup
% ==========================
function InitializeConditions(block)
    IC = block.DialogPrm(2).Data;

    % Convert angles and rates to radians
    P = IC.P*pi/180; Q = IC.Q*pi/180; R = IC.R*pi/180; 
    Phi = IC.Phi*pi/180; The = IC.The*pi/180; Psi = IC.Psi*pi/180;

    % Initial body velocities and position
    U = IC.U; V = IC.V; W = IC.W; 
    X = IC.X; Y = IC.Y; Z = IC.Z;

    init = [P, Q, R, Phi, The, Psi, U, V, W, X, Y, Z];

    % Set initial state values
    for i = 1:12
        block.OutputPort(i).Data = init(i);
        block.ContStates.Data(i) = init(i);
    end
end

% ==========================
% Output Update
% ==========================
function Outputs(block)
    for i = 1:12
        block.OutputPort(i).Data = block.ContStates.Data(i);
    end
end

% ==========================
% Continuous-Time Dynamics
% ==========================
function Derivatives(block)

    quad = block.DialogPrm(1).Data;

    % Retrieve state variables
    P = block.ContStates.Data(1);
    Q = block.ContStates.Data(2);
    R = block.ContStates.Data(3);
    Phi = block.ContStates.Data(4);
    The = block.ContStates.Data(5);
    Psi = block.ContStates.Data(6);
    U = block.ContStates.Data(7);
    V = block.ContStates.Data(8);
    W = block.ContStates.Data(9);
    X = block.ContStates.Data(10);
    Y = block.ContStates.Data(11);
    Z = block.ContStates.Data(12);

    % Motor speeds
    w1 = block.InputPort(1).Data;
    w2 = block.InputPort(2).Data;
    w3 = block.InputPort(3).Data;
    w4 = block.InputPort(4).Data;
    w  = [w1; w2; w3; w4];

    % External disturbance torques and forces
    Dist_tau = block.InputPort(5).Data(1:3);
    Dist_F   = block.InputPort(5).Data(4:6);

    % Gyroscopic moment due to rotating propellers
    tau_motorGyro = [Q*quad.Jm*2*pi/60*(-w1-w3+w2+w4);
                     P*quad.Jm*2*pi/60*(w1+w3-w2-w4);
                     0];

    % Total body moment and thrust force
    Mb = (quad.dctcq*(w.^2)) + tau_motorGyro + Dist_tau;
    Fb = [0; 0; sum(quad.ct*(w.^2))];

    % Angular velocity vector and skew-symmetric matrix
    omb_bi = [P; Q; R];
    OMb_bi = [ 0, -R,  Q;
               R,  0, -P;
              -Q,  P,  0];

    % Angular acceleration from Euler equation
    b_omdotb_bi = quad.Jbinv * (Mb - OMb_bi * quad.Jb * omb_bi);

    % Euler angle rate transformation matrix
    H_Phi = [1, tan(The)*sin(Phi), tan(The)*cos(Phi);
             0, cos(Phi),         -sin(Phi);
             0, sin(Phi)/cos(The), cos(Phi)/cos(The)];   

    % Euler angle derivatives
    Phidot = H_Phi * omb_bi;

    % Rotation matrix from body to inertial frame
    Rib = [cos(Psi)*cos(The), cos(Psi)*sin(The)*sin(Phi)-sin(Psi)*cos(Phi), cos(Psi)*sin(The)*cos(Phi)+sin(Psi)*sin(Phi);
           sin(Psi)*cos(The), sin(Psi)*sin(The)*sin(Phi)+cos(Psi)*cos(Phi), sin(Psi)*sin(The)*cos(Phi)-cos(Psi)*sin(Phi);
           -sin(The),         cos(The)*sin(Phi),                             cos(The)*cos(Phi)];

    Rbi = Rib';  % Inverse rotation (inertial to body)

    % Gravity and disturbance force in body frame
    ge = [0; 0; -quad.g];
    gb = Rbi * ge;
    Dist_Fb = Rbi * Dist_F;

    % Translational dynamics in body frame
    vb = [U; V; W];
    b_dv = (1/quad.mass)*Fb + gb + Dist_Fb - OMb_bi * vb;

    % Inertial position derivative
    i_dp = Rib * vb;

    % State derivatives
    dP = b_omdotb_bi(1);
    dQ = b_omdotb_bi(2);
    dR = b_omdotb_bi(3);
    dPhi   = Phidot(1);
    dTheta = Phidot(2);
    dPsi   = Phidot(3);
    dU = b_dv(1);
    dV = b_dv(2);
    dW = b_dv(3);
    dX = i_dp(1);
    dY = i_dp(2);
    dZ = i_dp(3);

    % Ground contact condition: prevent Z from going below 0
    if ((Z <= 0) && (dZ <= 0))
        dZ = 0;
        block.ContStates.Data(12) = 0;
    end

    % Return full derivative vector
    f = [dP dQ dR dPhi dTheta dPsi dU dV dW dX dY dZ]';
    block.Derivatives.Data = f;
end
