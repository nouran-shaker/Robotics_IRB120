function out = inverse_kinematics(T_desired, opts)
%INVERSE_KINEMATICS  Closed-form IK for the ABB IRB120 (all 8 solutions).
%
%   inverse_kinematics()              interactive: prompts for a target pose,
%                                     solves all 8 branches, lets you pick a
%                                     valid one to visualise.
%   inverse_kinematics(T)             solve for a given 4x4 pose T (mm).
%   out = inverse_kinematics(T, opts) headless/programmatic use.
%
%   opts fields (all optional):
%     .Plot        (true)  render workspace cloud + robot posture
%     .Interactive (true)  prompt when several valid solutions exist;
%                          if false, auto-pick the first valid one
%     .CloudPts    (2e5)   number of workspace-cloud sample points
%
%   The 8 solutions come from three independent binary choices:
%     SHOULDER : front reach (theta1) vs backswept (theta1 + pi)
%     ELBOW    : above vs below arm
%     WRIST    : wrist down vs wrist up
%   Each solution is verified by forward kinematics and screened against the
%   IRB120 joint limits. The chosen solution is exported to the base
%   workspace for the forward-kinematics / dynamics pipeline.

    if nargin < 1, clc; close all; end   % only clear when run interactively
    fprintf('====================================================\n');
    fprintf('   ABB IRB120 Inverse Kinematics Solver  (8 branches)\n');
    fprintf('====================================================\n\n');

    %% ---- options ----
    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'Plot'),        opts.Plot        = true;  end
    if ~isfield(opts,'Interactive'), opts.Interactive = true;  end
    if ~isfield(opts,'CloudPts'),    opts.CloudPts    = 2e5;   end

    %% ---- D-H parameters (mm) ----
    DH.d1 = 290;  DH.a2 = 270;  DH.a3 = 70;  DH.d4 = 302;  DH.d6 = 72;

    %% ---- joint limits (degrees) ----
    limits = [-165 165; -110 110; -110 70; -160 160; -120 120; -400 400];

    %% ---- target pose ----
    if nargin < 1 || isempty(T_desired)
        T_desired = get_user_input(DH);
    end
    if isempty(T_desired) || ~isequal(size(T_desired),[4 4])
        fprintf('Invalid input. Exiting...\n'); out = struct(); return;
    end

    fprintf('\n--- Desired Transformation Matrix ---\n');
    for row = 1:4
        fprintf('  [%10.4f %10.4f %10.4f %10.4f]\n', T_desired(row,:));
    end

    %% ---- the 8 configurations: columns = [SHOULDER ELBOW WRIST] ----
    config_list = [ 1 -1  1;    % Front  Above  Wrist-Down
                    1 -1 -1;    % Front  Above  Wrist-Up
                    1  1  1;    % Front  Below  Wrist-Down
                    1  1 -1;    % Front  Below  Wrist-Up
                   -1 -1  1;    % Back   Above  Wrist-Down
                   -1 -1 -1;    % Back   Above  Wrist-Up
                   -1  1  1;    % Back   Below  Wrist-Down
                   -1  1 -1];   % Back   Below  Wrist-Up
    config_names = {'Front Above - Wrist Down','Front Above - Wrist Up', ...
                    'Front Below - Wrist Down','Front Below - Wrist Up', ...
                    'Back Above - Wrist Down', 'Back Above - Wrist Up', ...
                    'Back Below - Wrist Down', 'Back Below - Wrist Up'};
    nC = size(config_list,1);

    %% ---- solve every branch ----
    all_theta_deg = nan(nC,6);
    all_theta_rad = nan(nC,6);
    pos_err   = nan(nC,1);
    rot_err   = nan(nC,1);
    is_valid  = false(nC,1);
    ik_solved = false(nC,1);
    fail_msg  = cell(nC,1);
    violated  = cell(nC,1);

    for k = 1:nC
        cfg.shoulder = config_list(k,1);
        cfg.elbow    = config_list(k,2);
        cfg.wrist    = config_list(k,3);
        try
            theta     = compute_ik(T_desired, cfg, DH);
            theta_deg = rad2deg(theta);
            ik_solved(k) = true;
            all_theta_rad(k,:) = theta';
            all_theta_deg(k,:) = theta_deg';

            T_check    = forward_kinematics(theta, DH);
            pos_err(k) = norm(T_check(1:3,4) - T_desired(1:3,4));
            rot_err(k) = norm(T_check(1:3,1:3) - T_desired(1:3,1:3), 'fro');

            bad = find(theta_deg(:) < limits(:,1) | theta_deg(:) > limits(:,2))';
            violated{k} = bad;
            is_valid(k) = isempty(bad);
        catch ME
            fail_msg{k} = ME.message;
        end
    end

    %% ---- display ----
    fprintf('\n====================================================\n');
    fprintf('   INVERSE KINEMATICS SOLUTIONS  (%d branches)\n', nC);
    fprintf('====================================================\n');
    for k = 1:nC
        fprintf('\n--- Configuration %d: %s ---\n', k, config_names{k});
        if ~ik_solved(k)
            fprintf('  FAILED: %s\n', fail_msg{k});
            continue;
        end
        td = all_theta_deg(k,:);
        fprintf('  Joint Angles (deg): [%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f]\n', td);
        fprintf('  Pos error: %.4e mm   Orient error: %.4e\n', pos_err(k), rot_err(k));
        if is_valid(k)
            fprintf('  >> VALID - all joints within limits\n');
        else
            fprintf('  >> INVALID - joints out of range:');
            for jj = violated{k}
                fprintf(' J%d(%.1f)', jj, td(jj));
            end
            fprintf('\n');
        end
    end

    valid_count = sum(is_valid);
    fprintf('\n====================================================\n');
    fprintf('  Valid solutions (within joint limits): %d / %d\n', valid_count, nC);
    fprintf('====================================================\n');

    %% ---- package output ----
    out = struct('config_names',{config_names}, 'config_list',config_list, ...
                 'theta_deg',all_theta_deg, 'theta_rad',all_theta_rad, ...
                 'pos_err',pos_err, 'rot_err',rot_err, 'is_valid',is_valid, ...
                 'T_desired',T_desired);

    if valid_count == 0
        fprintf('\n  No valid configuration found for this target.\n');
        return;
    end

    %% ---- choose a valid solution ----
    valid_idx = find(is_valid);
    nv = numel(valid_idx);
    if nv == 1 || ~opts.Interactive
        chosen = valid_idx(1);
        if opts.Interactive
            fprintf('\nOne valid solution found. Visualizing it...\n');
        end
    else
        fprintf('\nSelect which valid solution to visualize:\n');
        for i = 1:nv
            k = valid_idx(i);
            fprintf('  %d - %-26s [%.1f %.1f %.1f %.1f %.1f %.1f] deg\n', ...
                    i, config_names{k}, all_theta_deg(k,:));
        end
        sel = input(sprintf('\nEnter choice (1-%d): ', nv));
        if isempty(sel) || sel < 1 || sel > nv, sel = 1; end
        chosen = valid_idx(sel);
    end
    out.chosen = chosen;

    %% ---- export to base workspace (FK / dynamics pipeline) ----
    assignin('base', 'ik_q_deg',     all_theta_deg(chosen,:));
    assignin('base', 'ik_q_rad',     all_theta_rad(chosen,:));
    assignin('base', 'ik_config',    config_names{chosen});
    assignin('base', 'ik_T_desired', T_desired);

    %% ---- visualize ----
    if opts.Plot
        fprintf('\nGenerating workspace cloud (~%d pts)...\n', round(opts.CloudPts));
        cloud_pts = generate_workspace_cloud(DH, limits, round(opts.CloudPts));
        fprintf('Cloud ready. Plotting...\n');
        plot_solution(cloud_pts, all_theta_rad(chosen,:)', T_desired, ...
                      config_names{chosen}, all_theta_deg(chosen,:), DH);
    end

    fprintf('\n====================================================\n');
    fprintf('   IK Complete! Exported to workspace:\n');
    fprintf('   ik_q_deg / ik_q_rad / ik_config / ik_T_desired\n');
    fprintf('   Chosen: %s\n', config_names{chosen});
    fprintf('   >> Run forward_kinematics.m next\n');
    fprintf('====================================================\n');
end


%% ============================================================
%%  INVERSE KINEMATICS SOLVER  (geometric, Lee-Ziegler extended)
%% ============================================================
function theta = compute_ik(T, config, DH)
    d1=DH.d1; a2=DH.a2; a3=DH.a3; d4=DH.d4; d6=DH.d6;
    n=T(1:3,1); s=T(1:3,2); a=T(1:3,3); p=T(1:3,4);
    SHOULDER=config.shoulder; ELBOW=config.elbow; WRIST=config.wrist;

    % wrist centre (remove the tool offset along approach vector)
    p4=p-d6*a; px=p4(1); py=p4(2); pz=p4(3);

    % --- Joint 1: front reach vs backswept ---
    if SHOULDER > 0
        theta1 = atan2(py,px);
    else
        theta1 = atan2(py,px) + pi;
        theta1 = atan2(sin(theta1),cos(theta1));   % wrap to (-pi, pi]
    end

    % --- Joints 2 & 3: planar geometry (signed radius for backswept) ---
    r = SHOULDER*sqrt(px^2+py^2);            % signed planar radius
    R = sqrt(r^2+(pz-d1)^2);
    L_max=a2+sqrt(a3^2+d4^2); L_min=abs(a2-sqrt(a3^2+d4^2));
    if R>L_max || R<L_min
        error('Unreachable: R=%.1f mm (range [%.1f,%.1f])',R,L_min,L_max);
    end

    sin_alpha=(pz-d1)/R; cos_alpha=r/R;
    cos_beta=clamp_val((a2^2+R^2-a3^2-d4^2)/(2*a2*R));
    sin_beta=sqrt(1-cos_beta^2);
    sin_t2=-sin_alpha*cos_beta+ELBOW*cos_alpha*sin_beta;
    cos_t2= cos_alpha*cos_beta+ELBOW*sin_alpha*sin_beta;
    theta2=atan2(sin_t2,cos_t2);

    D=sqrt(a3^2+d4^2);
    cos_psi=clamp_val((a2^2+a3^2+d4^2-R^2)/(2*a2*D));
    sin_psi=sqrt(1-cos_psi^2);
    cos_phi=a3/D; sin_phi=d4/D;
    sin_t3= sin_phi*cos_psi-ELBOW*cos_phi*sin_psi;
    cos_t3=-cos_phi*cos_psi-ELBOW*sin_phi*sin_psi;
    theta3=atan2(sin_t3,cos_t3);

    % --- Wrist: orientation decoupling ---
    T1=dh(theta1,-pi/2,0,d1); T2=dh(theta2,0,a2,0); T3=dh(theta3,-pi/2,a3,0);
    T30=T1*T2*T3;
    z3=T30(1:3,3); x3=T30(1:3,1); y3=T30(1:3,2);

    z3xa=cross(z3,a); nc=norm(z3xa);
    if nc<1e-6
        theta4=0;                            % degenerate (a parallel to z3)
    else
        z4=z3xa/nc;
        sdz=dot(s,z4);
        if abs(sdz)>1e-6, Omega=sign(sdz); else, Omega=sign(dot(n,z4)); end
        M=WRIST*sign(Omega);
        theta4=atan2(-M*dot(z4,x3), M*dot(z4,y3));
    end

    T4=dh(theta4,-pi/2,0,d4); T40=T30*T4;
    x4=T40(1:3,1); y4=T40(1:3,2);
    theta5=atan2(dot(a,x4),-dot(a,y4));

    T5=dh(theta5,pi/2,0,0); T50=T40*T5;
    y5=T50(1:3,2);
    theta6=atan2(dot(n,y5),dot(s,y5));

    theta=[theta1;theta2;theta3;theta4;theta5;theta6];
end


%% ============================================================
%%  FORWARD KINEMATICS  (verification)
%% ============================================================
function T = forward_kinematics(theta, DH)
    T=dh(theta(1),-pi/2,0,DH.d1)*dh(theta(2),0,DH.a2,0)* ...
      dh(theta(3),-pi/2,DH.a3,0)*dh(theta(4),-pi/2,0,DH.d4)* ...
      dh(theta(5),pi/2,0,0)*dh(theta(6),0,0,DH.d6);
end

function A = dh(theta,alpha,a,d)
    ct=cos(theta);st=sin(theta);ca=cos(alpha);sa=sin(alpha);
    A=[ct,-ca*st,sa*st,a*ct; st,ca*ct,-sa*ct,a*st; 0,sa,ca,d; 0,0,0,1];
end

function v = clamp_val(v)
    v=max(-1,min(1,v));
end


%% ============================================================
%%  WORKSPACE CLOUD
%% ============================================================
function pts = generate_workspace_cloud(DH, limits_deg, N)
    if nargin<3, N=200000; end
    lim = deg2rad(limits_deg);
    q   = zeros(N,6);
    for j = 1:5                              % joint 6 does not move the tool tip
        q(:,j) = lim(j,1) + (lim(j,2)-lim(j,1)).*rand(N,1);
    end
    pts = zeros(N,3);
    for i = 1:N
        T = fk_cloud(q(i,:), DH);
        pts(i,:) = T(1:3,4)';
    end
end

function T = fk_cloud(q, DH)
    T = dh(q(1),-pi/2,0,DH.d1)*dh(q(2),0,DH.a2,0)* ...
        dh(q(3),-pi/2,DH.a3,0)*dh(q(4),-pi/2,0,DH.d4)* ...
        dh(q(5),pi/2,0,0)*dh(q(6),0,0,DH.d6);
end


%% ============================================================
%%  VISUALIZATION
%% ============================================================
function plot_solution(cloud_pts, theta_rad, T_desired, cfg_name, ang, DH)
    figure('Name','ABB IRB120 - IK Solution + Workspace', ...
           'Position',[50 50 900 750], 'Color','w');

    scatter3(cloud_pts(:,1), cloud_pts(:,2), cloud_pts(:,3), ...
             1, [0.60 0.92 0.98], '.', 'MarkerEdgeAlpha', 0.18);
    hold on;
    plot_robot(theta_rad, DH);

    p = T_desired(1:3,4);
    plot3(p(1),p(2),p(3),'g*','MarkerSize',14,'LineWidth',2,'DisplayName','Target');
    T_ee = forward_kinematics(theta_rad, DH);
    plot3(T_ee(1,4),T_ee(2,4),T_ee(3,4),'ro','MarkerSize',10, ...
          'MarkerFaceColor','r','LineWidth',1.5,'DisplayName','Tool Tip');

    legend({'Workspace Cloud','Robot Links','Target','Tool Tip'}, ...
           'Location','northeast','FontSize',10,'Color',[0.97 0.97 0.97]);

    angle_str = sprintf(['\\theta_1=%.1f  \\theta_2=%.1f  \\theta_3=%.1f  ' ...
                         '\\theta_4=%.1f  \\theta_5=%.1f  \\theta_6=%.1f'], ang);
    title({['ABB IRB120 - ' cfg_name], angle_str},'FontSize',12,'FontWeight','bold');
    xlabel('X (mm)','FontSize',11); ylabel('Y (mm)','FontSize',11); zlabel('Z (mm)','FontSize',11);
    grid on; axis equal; box on; view(135,25); set(gca,'FontSize',10);
end

function plot_robot(theta, DH)
    d1=DH.d1;a2=DH.a2;a3=DH.a3;d4=DH.d4;d6=DH.d6;
    T0=eye(4);
    T1=T0*dh(theta(1),-pi/2,0,d1); T2=T1*dh(theta(2),0,a2,0);
    T3=T2*dh(theta(3),-pi/2,a3,0); T4=T3*dh(theta(4),-pi/2,0,d4);
    T5=T4*dh(theta(5),pi/2,0,0);   T6=T5*dh(theta(6),0,0,d6);
    pts=[T0(1:3,4),T1(1:3,4),T2(1:3,4),T3(1:3,4),T4(1:3,4),T5(1:3,4),T6(1:3,4)];
    plot3(pts(1,:),pts(2,:),pts(3,:),'-','Color',[0.10 0.10 0.10],'LineWidth',4);
    plot3(pts(1,:),pts(2,:),pts(3,:),'o','MarkerSize',9, ...
          'MarkerFaceColor',[0.85 0.15 0.15],'MarkerEdgeColor','k','LineWidth',1);
    fill3([-80 80 80 -80],[-80 -80 80 80],[0 0 0 0], ...
          [0.55 0.55 0.55],'FaceAlpha',0.45,'EdgeColor','k');
    sc=55; draw_frame(T0,sc); draw_frame(T6,sc);
    max_r=d1+a2+sqrt(a3^2+d4^2)+d6;
    axis([-max_r max_r -max_r max_r -120 max_r*1.15]);
    set(gca,'FontSize',10);
end

function draw_frame(T,sc)
    o=T(1:3,4);
    quiver3(o(1),o(2),o(3),T(1,1)*sc,T(2,1)*sc,T(3,1)*sc,'r','LineWidth',2,'MaxHeadSize',0.6);
    quiver3(o(1),o(2),o(3),T(1,2)*sc,T(2,2)*sc,T(3,2)*sc,'g','LineWidth',2,'MaxHeadSize',0.6);
    quiver3(o(1),o(2),o(3),T(1,3)*sc,T(2,3)*sc,T(3,3)*sc,'b','LineWidth',2,'MaxHeadSize',0.6);
end


%% ============================================================
%%  USER INPUT
%% ============================================================
function T = get_user_input(DH) %#ok<INUSD>
    fprintf('Select input mode:\n');
    fprintf('  1 - Enter position (x, y, z) with default orientation\n');
    fprintf('  2 - Enter full 4x4 transformation matrix\n');
    fprintf('  3 - Enter position + Euler ZYZ orientation\n');
    fprintf('  4 - Use a pre-tested valid example\n');
    mode = input('\nEnter your choice (1-4): ');
    switch mode
        case 1
            fprintf('\nEnter desired end-effector position (mm):\n');
            x=input('  x = '); y=input('  y = '); z=input('  z = ');
            T=eye(4); T(1:3,4)=[x;y;z];
        case 2
            fprintf('\nEnter each row of the 4x4 matrix (space-separated):\n');
            T=zeros(4,4);
            for i=1:4
                row_str=input(sprintf('  Row %d: ',i),'s');
                T(i,:)=str2num(row_str); %#ok<ST2NM>
            end
        case 3
            fprintf('\nEnter desired end-effector position (mm):\n');
            x=input('  x = '); y=input('  y = '); z=input('  z = ');
            fprintf('\nEnter Euler ZYZ angles (degrees):\n');
            a1=deg2rad(input('  Alpha (rotation about Z): '));
            a2=deg2rad(input('  Beta  (rotation about Y): '));
            a3=deg2rad(input('  Gamma (rotation about Z): '));
            Rz1=[cos(a1) -sin(a1) 0;sin(a1) cos(a1) 0;0 0 1];
            Ry =[cos(a2) 0 sin(a2);0 1 0;-sin(a2) 0 cos(a2)];
            Rz2=[cos(a3) -sin(a3) 0;sin(a3) cos(a3) 0;0 0 1];
            T=[Rz1*Ry*Rz2,[x;y;z];0 0 0 1];
        case 4
            T=pick_example();
        otherwise
            T=[];
    end
end

function T = pick_example()
    fprintf('\nPre-tested valid examples:\n');
    fprintf('  1 - Front reach   (x=456, y=369, z=394)\n');
    fprintf('  2 - High reach    (x=-77, y=7,   z=900)\n');
    fprintf('  3 - Side reach    (x=-355,y=-174, z=374)\n');
    fprintf('  4 - Low side      (x=299, y=446,  z=42 )\n');
    ch = input('\nSelect example (1-4): ');
    switch ch
        case 1
            T=[ 0.1283, 0.1969, 0.9720, 455.98;
                0.8884, 0.4128,-0.2009, 368.52;
               -0.4408, 0.8893,-0.1219, 394.35;
                     0,      0,      0,      1];
        case 2
            T=[ 0.0251, 0.8557,-0.5169, -76.81;
               -0.6128, 0.4217, 0.6683,   7.11;
                0.7898, 0.3000, 0.5349, 899.74;
                     0,      0,      0,      1];
        case 3
            T=[-0.5594, 0.8036, 0.2032,-355.22;
                0.7906, 0.5909,-0.1605,-174.35;
               -0.2491, 0.0708,-0.9659, 373.49;
                     0,      0,      0,      1];
        case 4
            T=[ 0.0682, 0.6780,-0.7319, 299.04;
                0.5376, 0.5930, 0.5994, 445.68;
                0.8404,-0.4343,-0.3240,  41.87;
                     0,      0,      0,      1];
        otherwise
            fprintf('Invalid, using example 1.\n');
            T=[ 0.1283, 0.1969, 0.9720, 455.98;
                0.8884, 0.4128,-0.2009, 368.52;
               -0.4408, 0.8893,-0.1219, 394.35;
                     0,      0,      0,      1];
    end
end
