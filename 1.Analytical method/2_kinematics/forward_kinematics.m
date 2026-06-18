function forward_kinematics()
clc; close all;

fprintf('====================================================\n');
fprintf('   ABB IRB120 - Forward Kinematics\n');
fprintf('====================================================\n\n');

%% ============================================================
%  READ FROM IK WORKSPACE
% ============================================================
if ~evalin('base','exist(''ik_q_deg'',''var'')')
    error('ik_q_deg not found in workspace. Run inverse_kinematics first.');
end

q_deg = evalin('base', 'ik_q_deg');
q     = deg2rad(q_deg);

fprintf('Reading joint angles from IK workspace...\n');
if evalin('base','exist(''ik_config'',''var'')')
    fprintf('IK Configuration: %s\n', evalin('base','ik_config'));
end
fprintf('q_deg = ['); fprintf('%.3f  ', q_deg); fprintf(']\n\n');

%% ============================================================
%  PARAMETERS
% ============================================================
params = robot_params();

%% ============================================================
%  JOINT LIMIT CHECK
% ============================================================
check_joint_limits(q_deg, params);

%% ============================================================
%  FORWARD KINEMATICS
% ============================================================
[T_0_6, T_frames] = fk_compute(q, params);
position = T_0_6(1:3,4);

%% ============================================================
%  OUTPUT
% ============================================================
fprintf('====================================================\n');
fprintf('   RESULTS\n');
fprintf('====================================================\n');
fprintf('Transformation Matrix T_0_6:\n');
for i = 1:4
    fprintf('  ['); fprintf('%10.5f ', T_0_6(i,:)); fprintf(']\n');
end
fprintf('\nEnd-Effector Position:\n');
fprintf('  X = %.5f m\n', position(1));
fprintf('  Y = %.5f m\n', position(2));
fprintf('  Z = %.5f m\n', position(3));

%% ============================================================
%  PLOT
% ============================================================
plot_robot(T_frames, q_deg);

%% ============================================================
%  EXPORT TO WORKSPACE → Jacobian will read these
% ============================================================
assignin('base', 'q',        q);
assignin('base', 'q_deg',    q_deg);
assignin('base', 'T_frames', T_frames);
assignin('base', 'T_0_6',    T_0_6);

fprintf('\n====================================================\n');
fprintf('   FK Complete! Exported to workspace:\n');
fprintf('   q         - joint angles (radians) [1x6]\n');
fprintf('   q_deg     - joint angles (degrees) [1x6]\n');
fprintf('   T_frames  - all 7 frame transforms {7x1 cell}\n');
fprintf('   T_0_6     - end-effector transform [4x4]\n');
fprintf('   >> Run jacobian.m next\n');
fprintf('====================================================\n');

end


%% ============================================================
%%  ROBOT PARAMETERS
%% ============================================================
function params = robot_params()
    params.d1 = 0.290; params.a2 = 0.270;
    params.a3 = 0.070; params.d4 = 0.302; params.d6 = 0.072;
    params.joint_limits_deg = [-165,165; -110,110; -110,70;
                               -160,160; -120,120; -400,400];
    params.joint_limits = deg2rad(params.joint_limits_deg);
end


%% ============================================================
%%  JOINT LIMIT CHECK
%% ============================================================
function check_joint_limits(q_deg, params)
    fprintf('Joint Limit Validation:\n');
    failed = false;
    for j = 1:6
        lo = params.joint_limits_deg(j,1);
        hi = params.joint_limits_deg(j,2);
        v  = q_deg(j);
        if v < lo || v > hi
            fprintf('  [FAIL] Joint %d = %8.2f deg  limits [%.1f, %.1f]\n',j,v,lo,hi);
            failed = true;
        else
            fprintf('  [ OK ] Joint %d = %8.2f deg  limits [%.1f, %.1f]\n',j,v,lo,hi);
        end
    end
    if failed
        error('Joint limit violation detected. Fix IK output before running FK.');
    end
    fprintf('\n');
end


%% ============================================================
%%  FORWARD KINEMATICS
%% ============================================================
function [T_0_6, T_frames] = fk_compute(q, DH)
    T_frames    = cell(7,1);
    T_frames{1} = eye(4);
    T_frames{2} = T_frames{1} * dh_mat(q(1),-pi/2, 0,    DH.d1);
    T_frames{3} = T_frames{2} * dh_mat(q(2), 0,    DH.a2,0    );
    T_frames{4} = T_frames{3} * dh_mat(q(3),-pi/2, DH.a3,0    );
    T_frames{5} = T_frames{4} * dh_mat(q(4),-pi/2, 0,    DH.d4);
    T_frames{6} = T_frames{5} * dh_mat(q(5), pi/2, 0,    0    );
    T_frames{7} = T_frames{6} * dh_mat(q(6), 0,    0,    DH.d6);
    T_0_6 = T_frames{7};
end


%% ============================================================
%%  DH MATRIX
%% ============================================================
function A = dh_mat(theta, alpha, a, d)
    ct=cos(theta); st=sin(theta);
    ca=cos(alpha); sa=sin(alpha);
    A=[ct,-ca*st,sa*st,a*ct; st,ca*ct,-sa*ct,a*st; 0,sa,ca,d; 0,0,0,1];
end


%% ============================================================
%%  PLOT
%% ============================================================
function plot_robot(T_frames, q_deg)
    n   = length(T_frames);
    pts = zeros(3,n);
    for i = 1:n
        pts(:,i) = T_frames{i}(1:3,4);
    end

    figure('Name','ABB IRB120 - Forward Kinematics', ...
           'Color','w','Position',[100 100 800 650]);

    plot3(pts(1,:),pts(2,:),pts(3,:),'-o', ...
          'LineWidth',3,'MarkerSize',8, ...
          'Color',[0.15 0.35 0.80], ...
          'MarkerFaceColor',[0.90 0.20 0.20], ...
          'MarkerEdgeColor','k');

    grid on; axis equal; box on;
    xlabel('X (m)','FontSize',12); ylabel('Y (m)','FontSize',12);
    zlabel('Z (m)','FontSize',12);
    title(sprintf('FK — [%.1f  %.1f  %.1f  %.1f  %.1f  %.1f] deg', q_deg), ...
          'FontSize',13,'FontWeight','bold','Interpreter','none');
    view(45,30); set(gca,'FontSize',10);
end
