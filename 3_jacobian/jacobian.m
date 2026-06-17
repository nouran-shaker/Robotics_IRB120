function J = jacobian(q, T_frames, q_dot, opts)
%JACOBIAN  Geometric Jacobian, velocity mapping (both ways), singularity
%          classification, manipulability and static force duality for the
%          ABB IRB120.
%
%   jacobian()                          reads q, T_frames (and optional
%                                       q_dot) from the base workspace, as
%                                       produced by forward_kinematics.
%   J = jacobian(q, T_frames, q_dot)    run directly with given inputs.
%   J = jacobian(..., opts)             opts.Plots / opts.Verbose (logical),
%                                       opts.x_dot_des (6x1 desired EE twist).

clc;
if nargin < 4, opts = struct(); end
if ~isfield(opts,'Plots'),   opts.Plots   = true;  end
if ~isfield(opts,'Verbose'), opts.Verbose = true;  end

if opts.Verbose
    fprintf('====================================================\n');
    fprintf('   ABB IRB120 - Geometric Jacobian & Velocity Analysis\n');
    fprintf('====================================================\n\n');
end

%% ---- inputs (args or base workspace) ----
if nargin < 1 || isempty(q)
    if ~evalin('base','exist(''q'',''var'')') || ~evalin('base','exist(''T_frames'',''var'')')
        error('jacobian:noState','q or T_frames not found. Run forward_kinematics first.');
    end
    q        = evalin('base','q');
    T_frames = evalin('base','T_frames');
end
q = q(:)';

% joint velocities / accelerations (defaults preserved from original)
if nargin < 3 || isempty(q_dot)
    if nargin < 1 && evalin('base','exist(''q_dot'',''var'')')
        q_dot = evalin('base','q_dot');
    else
        q_dot = [0.2; 0.2; 0.0; 0.5; 0.2; 0.0];
    end
end
q_dot  = q_dot(:);
q_ddot = [0.1; 0.1; 0.0; 0.2; 0.1; 0.0];

% DH lengths (m) used only for the singularity geometry
DH.d1=0.290; DH.a2=0.270; DH.a3=0.070; DH.d4=0.302; DH.d6=0.072;

%% ---- origins and z-axes from FK frames ----
O = zeros(3,7); Z = zeros(3,7);
for i = 1:7
    O(:,i) = T_frames{i}(1:3,4);
    Z(:,i) = T_frames{i}(1:3,3);
end

%% ---- geometric Jacobian [6x6] ----
%  J(:,i) = [ z_{i-1} x (o_n - o_{i-1}) ;  z_{i-1} ]
J = zeros(6,6);
for i = 1:6
    J(:,i) = [cross(Z(:,i), O(:,7) - O(:,i)); Z(:,i)];
end

%% ---- forward velocity map  (joint -> task) ----
x_dot = J * q_dot;
v = x_dot(1:3);     % linear  (m/s)
w = x_dot(4:6);     % angular (rad/s)

%% ---- inverse velocity map  (task -> joint) ----  [Part C: <-> ]
% Recover joint rates from a desired end-effector twist. Uses a direct
% solve when well-conditioned, and damped least squares (DLS) near a
% singularity so the result stays bounded.
if isfield(opts,'x_dot_des'), x_des = opts.x_dot_des(:); else, x_des = [0;0;0.1;0;0;0]; end
lambda    = 0.01;                                % DLS damping factor
near_sing = min(svd(J)) < 1e-6;
if near_sing
    q_dot_direct = nan(6,1);                     % exact inverse undefined here
    res_fwd      = NaN;
else
    q_dot_direct = J \ x_des;                     % exact inverse map
    res_fwd      = norm(J*(J\x_dot) - x_dot);     % round-trip residual
end
q_dot_dls = J' * ((J*J' + lambda^2*eye(6)) \ x_des);   % robust near singularities

%% ---- static force / torque duality  tau = J' * F ----
F_ext = [0;0;-10;0;0;0];                         % example 10 N pull in -Z
tau_static = J' * F_ext;

%% ---- singularity analysis (numeric) ----
[~,S,V]        = svd(J);
sv             = diag(S);
sigma_min      = min(sv);
manipulability = prod(sv);                        % = |det J| for square J
detJ           = det(J);
yoshikawa      = sqrt(max(det(J*J'),0));          % sqrt(det(J J^T))
r              = rank(J);
cond_num       = max(sv) / (sigma_min + eps);
threshold      = 1e-4;

%% ---- singularity CLASSIFICATION (analytic geometry) ----  [Part C]
% wrist centre = intersection of axes 4,5,6 (origin of frame 4)
wc        = O(:,5);
shoulder  = O(:,2);
R_reach   = norm(wc - shoulder);
L_max     = DH.a2 + sqrt(DH.a3^2 + DH.d4^2);
L_min     = abs(DH.a2 - sqrt(DH.a3^2 + DH.d4^2));
tol_ang   = deg2rad(3);
tol_len   = 0.010;
sing_wrist    = abs(sin(q(5)))           < sin(tol_ang);   % axes 4 & 6 aligned
sing_elbow    = abs(R_reach - L_max)     < tol_len || abs(R_reach - L_min) < tol_len;
sing_shoulder = norm(wc(1:2))            < tol_len;        % wrist over base axis

%% ---- display ----
if opts.Verbose
    fprintf('Joint velocities  q_dot  = [%s] rad/s\n',  sprintf('%.2f ',q_dot));
    fprintf('Joint accels      q_ddot = [%s] rad/s^2\n\n', sprintf('%.2f ',q_ddot));

    fprintf('====================================================\n');
    fprintf('   GEOMETRIC JACOBIAN  [6x6]\n');
    fprintf('====================================================\n');
    rl = {'Vx','Vy','Vz','Wx','Wy','Wz'};
    for i = 1:6
        fprintf('  %s: %s\n', rl{i}, sprintf('%9.4f ', J(i,:)));
    end

    fprintf('\n====================================================\n');
    fprintf('   FORWARD MAP  (joint -> task):  x_dot = J q_dot\n');
    fprintf('====================================================\n');
    fprintf('  Linear  (m/s):   Vx=%8.4f  Vy=%8.4f  Vz=%8.4f\n', v);
    fprintf('  Angular (rad/s): Wx=%8.4f  Wy=%8.4f  Wz=%8.4f\n', w);

    fprintf('\n====================================================\n');
    fprintf('   INVERSE MAP  (task -> joint):  q_dot = J^{-1} x_dot\n');
    fprintf('====================================================\n');
    fprintf('  Desired EE twist: [%s]\n', sprintf('%.3f ', x_des));
    if near_sing
        fprintf('  q_dot (direct) : n/a  (Jacobian singular - exact inverse undefined)\n');
    else
        fprintf('  q_dot (direct) : [%s] rad/s\n', sprintf('%8.4f ', q_dot_direct));
        fprintf('  Round-trip residual ||J(J\\x)-x|| = %.3e\n', res_fwd);
    end
    fprintf('  q_dot (DLS)    : [%s] rad/s   (lambda=%.3f)\n', sprintf('%8.4f ', q_dot_dls), lambda);

    fprintf('\n====================================================\n');
    fprintf('   STATIC FORCE DUALITY:  tau = J^T F\n');
    fprintf('====================================================\n');
    fprintf('  EE wrench F   : [%s]\n', sprintf('%.1f ', F_ext));
    fprintf('  Joint torques : [%s] N.m\n', sprintf('%8.4f ', tau_static));

    fprintf('\n====================================================\n');
    fprintf('   SINGULARITY ANALYSIS\n');
    fprintf('====================================================\n');
    fprintf('  Rank:                 %d\n', r);
    fprintf('  Singular values:      %s\n', sprintf('%.4f ', sv'));
    fprintf('  Min singular value:   %.6f\n', sigma_min);
    fprintf('  det(J):               %.6f\n', detJ);
    fprintf('  Manipulability |det|: %.6f\n', manipulability);
    fprintf('  Yoshikawa sqrt(JJ''): %.6f\n', yoshikawa);
    if sigma_min > 1e-10
        fprintf('  Condition number:     %.2f\n', cond_num);
    else
        fprintf('  Condition number:     Inf (singular)\n');
    end

    if sigma_min < threshold
        nv = V(:,end);
        [~,idx] = max(abs(nv));
        fprintf('\n  >> SINGULARITY DETECTED\n');
        fprintf('  Null-space direction: [%s]\n', sprintf('%.4f ', nv'));
        fprintf('  Most affected joint:  Joint %d\n', idx);
    else
        fprintf('\n  >> Non-singular configuration. All 6 DOF available.\n');
    end
    fprintf('  --- type check (geometry) ---\n');
    fprintf('  Wrist  (theta5 ~ 0/180): %s\n', yn(sing_wrist));
    fprintf('  Elbow  (R = %.3f m near [%.3f,%.3f]): %s\n', R_reach, L_min, L_max, yn(sing_elbow));
    fprintf('  Shoulder (wrist over base axis): %s\n', yn(sing_shoulder));
    fprintf('====================================================\n');
end

%% ---- visualization ----
if opts.Plots
    plot_jacobian(O, Z, v, w, r, sigma_min, J);
end

%% ---- export ----
assignin('base','J',J);
assignin('base','q_dot',q_dot);
assignin('base','q_ddot',q_ddot);
assignin('base','x_dot',x_dot);
assignin('base','manipulability',manipulability);

if opts.Verbose
    fprintf('\n====================================================\n');
    fprintf('   Jacobian Complete! Exported: J, q_dot, q_ddot, x_dot, manipulability\n');
    fprintf('   >> Run irb120_dynamics.m next\n');
    fprintf('====================================================\n');
end
end


%% ============================================================
function s = yn(tf)
    if tf, s = 'YES'; else, s = 'no'; end
end


%% ============================================================
%%  VISUALIZATION  (links, velocity vectors, axes, manipulability ellipsoid)
%% ============================================================
function plot_jacobian(O, Z, v, w, r, sigma_min, J)
    figure('Name','Jacobian Velocity Analysis','Position',[100 100 900 700],'Color','w');
    hold on; grid on; axis equal; box on;

    plot3(O(1,:),O(2,:),O(3,:),'-o','LineWidth',2.5,'MarkerSize',8, ...
          'Color',[0.15 0.35 0.80],'MarkerFaceColor',[0.90 0.20 0.20],'MarkerEdgeColor','k');

    scale_v = 0.3; scale_w = 0.1;
    ee = O(:,7);
    quiver3(ee(1),ee(2),ee(3), v(1)*scale_v,v(2)*scale_v,v(3)*scale_v, ...
            0,'r','LineWidth',3,'MaxHeadSize',0.5,'DisplayName','Linear Vel');
    quiver3(ee(1),ee(2),ee(3), w(1)*scale_w,w(2)*scale_w,w(3)*scale_w, ...
            0,'b','LineWidth',3,'MaxHeadSize',0.5,'DisplayName','Angular Vel');
    for i = 1:6
        quiver3(O(1,i),O(2,i),O(3,i), Z(1,i)*0.05,Z(2,i)*0.05,Z(3,i)*0.05, ...
                0,'g','LineWidth',1.5,'MaxHeadSize',0.5);
    end

    % translational velocity (manipulability) ellipsoid at the end-effector
    Jv = J(1:3,:); A = Jv*Jv';
    [Vec,Dg] = eig((A+A')/2);
    radii = sqrt(max(diag(Dg),0)) * scale_v;
    [xe,ye,ze] = ellipsoid(0,0,0,radii(1),radii(2),radii(3),18);
    P = Vec * [xe(:)'; ye(:)'; ze(:)'];
    Xe = reshape(P(1,:),size(xe)) + ee(1);
    Ye = reshape(P(2,:),size(ye)) + ee(2);
    Ze = reshape(P(3,:),size(ze)) + ee(3);
    surf(Xe,Ye,Ze,'FaceColor',[0.55 0.55 0.90],'FaceAlpha',0.18, ...
         'EdgeColor','none','DisplayName','Velocity Ellipsoid');

    legend({'Robot Links','Linear Velocity','Angular Velocity','Joint Axes','Velocity Ellipsoid'}, ...
           'Location','best','FontSize',10,'Interpreter','none');
    xlabel('X (m)','FontSize',12,'FontWeight','bold');
    ylabel('Y (m)','FontSize',12,'FontWeight','bold');
    zlabel('Z (m)','FontSize',12,'FontWeight','bold');
    title(sprintf('Jacobian Analysis   Rank=%d   sigma_min=%.4f', r, sigma_min), ...
          'FontSize',13,'FontWeight','bold','Interpreter','none');
    view(45,30); set(gca,'FontSize',11);
end
