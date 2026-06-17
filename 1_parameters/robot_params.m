function params = robot_params()
    % ====================================================================
    % ABB IRB 120-3/0.6 ROBOT PARAMETERS
    % ====================================================================
    % STANDARD DH CONVENTION (Classic Denavit-Hartenberg)
    % 
    % Reference: Product Specification ROBO149EN Rev. J (November 2019)
    % ====================================================================
    
    %% ====================================================================
    %  GEOMETRIC PARAMETERS (in METERS)
    %  ====================================================================
    
    % Link dimensions (convert mm to meters)
    params.d1 = 0.290;      % 290 mm - Base height
    params.a2 = 0.270;      % 270 mm - Upper arm length
    params.a3 = 0.070;      % 70 mm  - Forearm offset
    params.d4 = 0.302;      % 302 mm - Wrist length
    params.d6 = 0.072;      % 72 mm  - Tool flange distance
    
    %% ====================================================================
    %  DENAVIT-HARTENBERG (DH) PARAMETERS
    %  ====================================================================
    %  STANDARD DH CONVENTION
    %  ====================================================================
    
   params.DH = [
    % theta_offset   d              a              alpha
    0,               params.d1,     0,             -pi/2;   % Joint 1
    0,               0,             params.a2,      0;      % Joint 2  
    0,               0,             params.a3,     -pi/2;   % Joint 3
    0,               params.d4,     0,             -pi/2;   % Joint 4  
    0,               0,             0,              pi/2;   % Joint 5  
    0,               params.d6,     0,              0;      % Joint 6
];
    %% ====================================================================
    %  JOINT LIMITS (from datasheet)
    %  ====================================================================
    
    params.joint_limits_deg = [
        -165,   +165;   % Joint 1: ±165°
        -110,   +110;   % Joint 2: -110° to +110°
        -110,   +70;    % Joint 3: -110° to +70° (ASYMMETRIC)
        -160,   +160;   % Joint 4: ±160°
        -120,   +120;   % Joint 5: ±120°
        -400,   +400;   % Joint 6: ±400° (multi-turn capability)
    ];
    
    % Convert to radians
    params.joint_limits = params.joint_limits_deg * pi/180;
    
    %% ====================================================================
    %  MAXIMUM JOINT VELOCITIES (from datasheet)
    %  ====================================================================
    
    params.max_joint_vel_deg = [
        250;    % Joint 1: 250°/s
        250;    % Joint 2: 250°/s
        250;    % Joint 3: 250°/s
        320;    % Joint 4: 320°/s
        320;    % Joint 5: 320°/s
        420;    % Joint 6: 420°/s
    ];
    
    % Convert to radians/second
    params.max_joint_vel = params.max_joint_vel_deg * pi/180;
    
    %% ====================================================================
    %  PERFORMANCE SPECIFICATIONS (from datasheet)
    %  ====================================================================
    
    params.max_payload = 3.0;                % 3 kg (4 kg with vertical wrist)
    params.max_reach = 0.580;                % 580 mm horizontal reach
    params.position_repeatability = 0.01e-3; % 0.01 mm (±0.01 mm)
    params.robot_mass = 25;                  % 25 kg total robot weight
    
    % Expected TCP at zero position
    params.tcp_zero_position = [0.374; 0.000; 0.630];  % [X; Y; Z] in meters
    
    %% ====================================================================
    %  DYNAMIC PARAMETERS
    %  ====================================================================
    %  Link masses (estimated distribution of 25 kg total)
    
    params.masses = [
        7.0;    % Link 1: Base + Joint 1 motor (~28%)
        6.5;    % Link 2: Shoulder + Upper arm (~26%)
        4.0;    % Link 3: Elbow + Forearm (~16%)
        3.5;    % Link 4: Wrist assembly 1 (~14%)
        2.5;    % Link 5: Wrist assembly 2 (~10%)
        1.5;    % Link 6: Flange + Joint 6 motor (~6%)
    ];
    
    % Verify total mass
    assert(abs(sum(params.masses) - params.robot_mass) < 0.1, ...
        'Link masses do not sum to total robot mass!');
    
    %% Center of Mass Positions
    %  Position [x, y, z] of each link's COM relative to its joint frame (meters)
    
    params.com = [
        0.000,  0.000,  0.145;   % Link 1: Middle of base height
        0.135,  0.000,  0.000;   % Link 2: Middle of upper arm
        0.035,  0.000,  0.000;   % Link 3: Along forearm offset
        0.000,  0.000,  0.151;   % Link 4: Middle of wrist length
        0.000,  0.000,  0.000;   % Link 5: At joint center
        0.000,  0.000,  0.036;   % Link 6: Middle of flange
    ];
    
    %% Inertia Tensors
    %  3×3 inertia matrix for each link about its center of mass (kg·m²)
    %  Using simplified cylinder/box approximations
    
    % Link 1: Cylindrical base (vertical cylinder)
    r1 = 0.090;  % Radius
    h1 = params.d1;  % Height
    m1 = params.masses(1);
    params.I{1} = diag([
        (1/12)*m1*(3*r1^2 + h1^2);  % Ixx
        (1/12)*m1*(3*r1^2 + h1^2);  % Iyy
        (1/2)*m1*r1^2;              % Izz
    ]);
    
    % Link 2: Upper arm (horizontal cylinder)
    r2 = 0.050;  % Radius
    L2 = params.a2;  % Length
    m2 = params.masses(2);
    params.I{2} = diag([
        (1/12)*m2*(3*r2^2 + L2^2);  % Ixx
        (1/12)*m2*(3*r2^2 + L2^2);  % Iyy
        (1/2)*m2*r2^2;              % Izz
    ]);
    
    % Link 3: Forearm (short cylinder)
    r3 = 0.040;  % Radius
    L3 = params.a3;  % Length
    m3 = params.masses(3);
    params.I{3} = diag([
        (1/12)*m3*(3*r3^2 + L3^2);  % Ixx
        (1/12)*m3*(3*r3^2 + L3^2);  % Iyy
        (1/2)*m3*r3^2;              % Izz
    ]);
    
    % Link 4: Wrist segment 1 (vertical cylinder)
    r4 = 0.030;  % Radius
    h4 = params.d4;  % Height
    m4 = params.masses(4);
    params.I{4} = diag([
        (1/12)*m4*(3*r4^2 + h4^2);  % Ixx
        (1/12)*m4*(3*r4^2 + h4^2);  % Iyy
        (1/2)*m4*r4^2;              % Izz
    ]);
    
    % Link 5: Wrist segment 2 (sphere approximation)
    r5 = 0.025;  % Radius
    m5 = params.masses(5);
    params.I{5} = diag([
        (2/5)*m5*r5^2;  % Ixx
        (2/5)*m5*r5^2;  % Iyy
        (2/5)*m5*r5^2;  % Izz
    ]);
    
    % Link 6: Tool flange (short cylinder)
    r6 = 0.040;  % Radius
    h6 = params.d6;  % Height
    m6 = params.masses(6);
    params.I{6} = diag([
        (1/12)*m6*(3*r6^2 + h6^2);  % Ixx
        (1/12)*m6*(3*r6^2 + h6^2);  % Iyy
        (1/2)*m6*r6^2;              % Izz
    ]);
    
    %% ====================================================================
    %  GRAVITY
    %  ====================================================================
    
    params.g = 9.81;  % Gravitational acceleration (m/s²)
    params.g_vector = [0; 0; -params.g];  % Gravity vector (downward)
    
    
    %% ====================================================================
    %  ROBOT IDENTIFICATION
    %  ====================================================================
    
    params.model = 'IRB 120-3/0.6';
    params.manufacturer = 'ABB';
    params.num_axes = 6;
    params.DH_convention = 'Standard DH (Classic)';
    params.datasheet_ref = 'ROBO149EN Rev. J (November 2019)';
    
end