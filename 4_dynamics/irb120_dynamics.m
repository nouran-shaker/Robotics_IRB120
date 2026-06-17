function results = irb120_dynamics(q, q_dot, q_ddot, varargin)
%IRB120_DYNAMICS  Rigid-body inverse-dynamic model of the ABB IRB120.
%
%   results = IRB120_DYNAMICS(q, q_dot, q_ddot) evaluates the joint-space
%   equation of motion for the 6-DOF ABB IRB120 manipulator:
%
%       tau = M(q)*q_ddot + C(q,q_dot)*q_dot + G(q)
%
%   M(q)   - configuration-dependent inertia matrix      [6x6] kg.m^2
%   C(.,.) - Coriolis/centrifugal matrix (Christoffel)   [6x6]
%   G(q)   - gravity torque vector                        [6x1] N.m
%   tau    - required joint torques                       [6x1] N.m
%
%   Inputs are joint angles in RADIANS, velocities in rad/s and
%   accelerations in rad/s^2 (row or column vectors, length 6). The
%   kinematic conventions (standard DH, frame ordering) are identical to
%   the companion Jacobian module so the two stay consistent.
%
%   Name-value options:
%     'Plots'   (true)   render the four diagnostic figures
%     'Verbose' (true)   print formatted tables to the command window
%     'Export'  (true)   push tau / M / C / G to the base workspace
%     'Step'    (1e-5)   central-difference step for Christoffel symbols
%
%   With no input arguments the state is read from the base workspace
%   (q, q_dot, q_ddot) for backward compatibility with the Jacobian script.
%
%   The returned struct contains M, C, G, tau, Mdot and a 'checks'
%   sub-struct holding the verification metrics.
%
%   Example:
%       q = deg2rad([10 -30 40 20 30 -15]);
%       r = irb120_dynamics(q, zeros(1,6), zeros(1,6));   % static hold

% ---- state acquisition -------------------------------------------------
if nargin < 3
    need = {'q','q_dot','q_ddot'};
    for k = 1:numel(need)
        if ~evalin('base', sprintf('exist(''%s'',''var'')', need{k}))
            error('irb120_dynamics:noState', ...
                  '%s not found in base workspace. Run the Jacobian script first.', need{k});
        end
    end
    q      = evalin('base','q');
    q_dot  = evalin('base','q_dot');
    q_ddot = evalin('base','q_ddot');
end

% ---- option parsing ----------------------------------------------------
opt = parse_options(varargin);

% ---- input validation --------------------------------------------------
q      = check_vec(q,      'q');
q_dot  = check_vec(q_dot,  'q_dot');
q_ddot = check_vec(q_ddot, 'q_ddot');
q_deg  = rad2deg(q);

p = get_params();
check_limits(q_deg, p);

if opt.Verbose
    banner('ABB IRB120  -  Inverse Dynamics    tau = M*q_ddot + C*q_dot + G');
    fprintf('  q      = [%s] deg\n',   join_num(q_deg));
    fprintf('  q_dot  = [%s] rad/s\n', join_num(q_dot));
    fprintf('  q_ddot = [%s] rad/s^2\n\n', join_num(q_ddot));
end

% ---- core dynamics -----------------------------------------------------
T_frames = compute_fk(q, p);                       % {1}=base ... {7}=EE
[Jv, Jw] = compute_link_jacobians(T_frames, p);
M_mat    = assemble_inertia(Jv, Jw, T_frames, p);

dMdq     = inertia_gradient(q, p, opt.Step);        % central differences
C_mat    = coriolis_from_gradient(dMdq, q_dot);     % Christoffel symbols
M_dot    = weighted_sum(dMdq, q_dot);               % d/dt M  along q_dot
G_vec    = compute_gravity(Jv, p);

tau = M_mat*q_ddot(:) + C_mat*q_dot(:) + G_vec;

% ---- verification ------------------------------------------------------
chk = run_checks(M_mat, C_mat, G_vec, tau, M_dot);
if opt.Verbose
    display_results(M_mat, C_mat, G_vec, tau, q_dot, q_ddot);
    display_checks(chk);
end

% ---- package + export --------------------------------------------------
results = struct('M', M_mat, 'C', C_mat, 'G', G_vec, 'tau', tau, ...
                 'Mdot', M_dot, 'q', q, 'q_dot', q_dot, 'q_ddot', q_ddot, ...
                 'checks', chk);

if opt.Export
    assignin('base','tau',   tau);
    assignin('base','M_mat', M_mat);
    assignin('base','C_mat', C_mat);
    assignin('base','G_vec', G_vec);
end

% ---- plots -------------------------------------------------------------
if opt.Plots
    th = theme();
    plot_component_breakdown(M_mat, C_mat, G_vec, q_ddot, q_dot, th);
    plot_inertia_heatmap(M_mat, th);
    plot_robot_posture(T_frames, q_deg, th);
end

if opt.Verbose
    banner('Dynamics complete.  Results returned as struct (and exported).');
end
end


%% ======================================================================
%   STATE / OPTION HELPERS
%  ======================================================================
function opt = parse_options(args)
    opt = struct('Plots', true, 'Verbose', true, 'Export', true, 'Step', 1e-5);
    if mod(numel(args), 2) ~= 0
        error('irb120_dynamics:badOptions','Options must be name-value pairs.');
    end
    for k = 1:2:numel(args)
        name = args{k};
        if ~isfield(opt, name)
            error('irb120_dynamics:unknownOption','Unknown option "%s".', name);
        end
        opt.(name) = args{k+1};
    end
end

function v = check_vec(v, name)
    if ~isnumeric(v) || numel(v) ~= 6 || ~isreal(v) || any(~isfinite(v(:)))
        error('irb120_dynamics:badInput', ...
              '%s must be a real, finite, 6-element vector.', name);
    end
    v = v(:)';                 % normalise to a row vector
end


%% ======================================================================
%   ROBOT PARAMETERS  (kinematics, masses, COMs, link inertias)
%  ======================================================================
function p = get_params()
    % --- DH lengths (m) ---
    p.d1 = 0.290; p.a2 = 0.270; p.a3 = 0.070; p.d4 = 0.302; p.d6 = 0.072;

    % --- joint range of motion (deg) ---
    p.limits_deg = [-165 165; -110 110; -110 70; -160 160; -120 120; -400 400];

    % --- gravity, masses (kg), COM offsets in each link frame (m) ---
    p.g_vec = [0; 0; -9.81];
    p.m     = [7.0; 6.5; 4.0; 3.5; 2.5; 1.5];
    p.com   = [0.000 0.000 0.145;
               0.135 0.000 0.000;
               0.035 0.000 0.000;
               0.000 0.000 0.151;
               0.000 0.000 0.000;
               0.000 0.000 0.036];

    % --- link inertia tensors about each COM (principal, link frame) ---
    % Links modelled as solid cylinders / rods; link 5 as a sphere (wrist).
    rad = [0.090 0.050 0.040 0.030 0.025 0.040];      % effective radii
    len = [p.d1  p.a2  p.a3  p.d4  0      p.d6 ];      % effective lengths
    p.I = cell(1,6);
    for i = 1:6
        m = p.m(i); r = rad(i); L = len(i);
        if i == 5                                     % spherical wrist link
            Ixx = (2/5)*m*r^2; Iyy = Ixx; Izz = Ixx;
        else                                          % cylinder about its axis
            Ixx = (1/12)*m*(3*r^2 + L^2);
            Iyy = Ixx;
            Izz = (1/2)*m*r^2;
        end
        p.I{i} = diag([Ixx Iyy Izz]);
    end
end


%% ======================================================================
%   FORWARD KINEMATICS  (standard DH, frame{1}=base, frame{7}=end-effector)
%  ======================================================================
function T = compute_fk(q, p)
    T    = cell(7,1);
    T{1} = eye(4);
    T{2} = T{1} * dh_mat(q(1), -pi/2, 0,    p.d1);
    T{3} = T{2} * dh_mat(q(2),  0,    p.a2, 0   );
    T{4} = T{3} * dh_mat(q(3), -pi/2, p.a3, 0   );
    T{5} = T{4} * dh_mat(q(4), -pi/2, 0,    p.d4);
    T{6} = T{5} * dh_mat(q(5),  pi/2, 0,    0   );
    T{7} = T{6} * dh_mat(q(6),  0,    0,    p.d6);
end

function A = dh_mat(theta, alpha, a, d)
    ct = cos(theta); st = sin(theta);
    ca = cos(alpha); sa = sin(alpha);
    A = [ct, -ca*st,  sa*st, a*ct;
         st,  ca*ct, -sa*ct, a*st;
          0,  sa,     ca,    d;
          0,  0,      0,     1];
end


%% ======================================================================
%   PER-LINK GEOMETRIC JACOBIANS (evaluated at each link COM)
%  ======================================================================
function [Jv, Jw] = compute_link_jacobians(T, p)
    n = 6; Jv = cell(n,1); Jw = cell(n,1);
    for i = 1:n
        r_com = T{i+1}(1:3,4) + T{i+1}(1:3,1:3)*p.com(i,:)';
        Jvi = zeros(3,n); Jwi = zeros(3,n);
        for j = 1:i
            z_j = T{j}(1:3,3);              % axis of joint j (z_{j-1})
            o_j = T{j}(1:3,4);
            Jvi(:,j) = cross(z_j, r_com - o_j);
            Jwi(:,j) = z_j;
        end
        Jv{i} = Jvi; Jw{i} = Jwi;
    end
end


%% ======================================================================
%   INERTIA MATRIX  M(q) = sum_i [ m_i Jv_i'Jv_i + Jw_i' R_i I_i R_i' Jw_i ]
%  ======================================================================
function M = assemble_inertia(Jv, Jw, T, p)
    M = zeros(6,6);
    for i = 1:6
        Ri = T{i+1}(1:3,1:3);
        Iw = Ri * p.I{i} * Ri';            % link inertia in world frame
        M  = M + p.m(i)*(Jv{i}'*Jv{i}) + Jw{i}'*Iw*Jw{i};
    end
    M = 0.5*(M + M');                       % enforce numerical symmetry
end

% inertia matrix evaluated directly from a configuration (for differencing)
function M = inertia_at(q, p)
    T        = compute_fk(q, p);
    [Jv, Jw] = compute_link_jacobians(T, p);
    M        = assemble_inertia(Jv, Jw, T, p);
end


%% ======================================================================
%   INERTIA GRADIENT  dM/dq_k  via CENTRAL differences  (O(h^2) accurate)
%  ======================================================================
function dMdq = inertia_gradient(q, p, h)
    n = 6; dMdq = zeros(n,n,n);
    for k = 1:n
        qp = q; qp(k) = qp(k) + h;
        qm = q; qm(k) = qm(k) - h;
        dMdq(:,:,k) = (inertia_at(qp,p) - inertia_at(qm,p)) / (2*h);
    end
end


%% ======================================================================
%   CORIOLIS MATRIX  C(i,j) = sum_k 1/2 (dM_ij/dq_k + dM_ik/dq_j - dM_jk/dq_i) q_dot_k
%  ======================================================================
function C = coriolis_from_gradient(dMdq, q_dot)
    n = size(dMdq,1); C = zeros(n);
    for i = 1:n
        for j = 1:n
            s = 0;
            for k = 1:n
                s = s + 0.5*(dMdq(i,j,k) + dMdq(i,k,j) - dMdq(j,k,i))*q_dot(k);
            end
            C(i,j) = s;
        end
    end
end

% d/dt M = sum_k (dM/dq_k) q_dot_k   (used for the passivity check)
function S = weighted_sum(dMdq, w)
    n = size(dMdq,1); S = zeros(n);
    for k = 1:n
        S = S + dMdq(:,:,k)*w(k);
    end
end


%% ======================================================================
%   GRAVITY VECTOR  G_i = -sum_{k>=i} m_k g' Jv_k(:,i)
%  ======================================================================
function G = compute_gravity(Jv, p)
    G = zeros(6,1);
    for i = 1:6
        for k = i:6
            G(i) = G(i) - p.m(k)*(p.g_vec' * Jv{k}(:,i));
        end
    end
end


%% ======================================================================
%   VERIFICATION
%  ======================================================================
function chk = run_checks(M, C, G, tau, M_dot)
    [gmax, gidx]   = max(abs(G));
    skew           = M_dot - 2*C;                    % must be skew-symmetric
    chk.sym_err    = norm(M - M', 'fro');
    chk.min_eig    = min(eig(M));
    chk.G1         = G(1);
    chk.grav_joint = gidx;
    chk.grav_max   = gmax;
    chk.passivity  = norm(skew + skew', 'fro');      % ||N + N'||  -> 0
    chk.static_res = norm(tau - G);                  % 0 only when q_dot=q_ddot=0
end


%% ======================================================================
%   CONSOLE OUTPUT
%  ======================================================================
function display_results(M, C, G, tau, q_dot, q_ddot)
    print_matrix('INERTIA MATRIX  M(q)  [6x6]  (kg.m^2)', M);
    print_matrix('CORIOLIS MATRIX  C(q,q_dot)  [6x6]',   C);

    banner('GRAVITY  G(q)   /   TORQUE  tau   (N.m)');
    Mc = M*q_ddot(:); Cc = C*q_dot(:);
    fprintf('  %-8s %12s %12s %12s %12s\n','Joint','M*q_ddot','C*q_dot','G','tau');
    fprintf('  %s\n', repmat('-',1,60));
    for i = 1:6
        fprintf('  J%-7d %12.5f %12.5f %12.5f %12.5f\n', i, Mc(i), Cc(i), G(i), tau(i));
    end
    fprintf('\n');
end

function display_checks(c)
    banner('VERIFICATION CHECKS');
    fprintf('  1) M symmetry error        = %.2e   (expect < 1e-10)\n', c.sym_err);
    fprintf('  2) M min eigenvalue        = %.6f   (must be > 0)\n',     c.min_eig);
    fprintf('  3) G(1) about vertical J1  = %.2e   (expect ~ 0)\n',      c.G1);
    fprintf('  4) Largest gravity torque  = J%d  (%.4f N.m)\n',          c.grav_joint, c.grav_max);
    fprintf('  5) Passivity ||Mdot-2C+()''|| = %.2e   (expect < 1e-6)\n', c.passivity);
    fprintf('  %s\n', repmat('-',1,52));
    pass = c.sym_err < 1e-10 && c.min_eig > 0 && c.passivity < 1e-6;
    if pass
        fprintf('  STATUS: ALL PHYSICAL CHECKS PASSED\n');
    else
        fprintf('  STATUS: CHECK FAILURE - review model\n');
    end
    fprintf('\n');
end

function print_matrix(ttl, A)
    banner(ttl);
    for i = 1:size(A,1)
        fprintf('  '); fprintf('%10.5f ', A(i,:)); fprintf('\n');
    end
    fprintf('\n');
end

function banner(txt)
    bar = repmat('=',1,68);
    fprintf('%s\n   %s\n%s\n', bar, txt, bar);
end

function s = join_num(v)
    s = sprintf('%8.3f', v(1));
    for k = 2:numel(v); s = [s sprintf('  %8.3f', v(k))]; end %#ok<AGROW>
end


%% ======================================================================
%   JOINT-LIMIT CHECK
%  ======================================================================
function check_limits(q_deg, p)
    bad = q_deg < p.limits_deg(:,1)' | q_deg > p.limits_deg(:,2)';
    if any(bad)
        j = find(bad,1);
        error('irb120_dynamics:limit', ...
              'Joint %d = %.2f deg violates limits [%.1f, %.1f].', ...
              j, q_deg(j), p.limits_deg(j,1), p.limits_deg(j,2));
    end
end


%% ======================================================================
%   PLOTTING THEME + FIGURES
%  ======================================================================
function th = theme()
    th.blue = [0.20 0.47 0.75];
    th.red  = [0.80 0.22 0.22];
    th.grn  = [0.10 0.68 0.38];
    th.gold = [1.00 0.80 0.00];
    th.ink  = [0.12 0.12 0.12];
    th.font = 11;
end

function plot_component_breakdown(M, C, G, q_ddot, q_dot, th)
    figure('Name','Torque Components','Position',[800 400 720 430],'Color','w');
    ax = axes; hold(ax,'on'); x = 1:6; w = 0.24;
    b1 = bar(ax,x-w, M*q_ddot(:), w,'FaceColor',th.blue,'EdgeColor','k','DisplayName','M\ddotq (inertial)');
    b2 = bar(ax,x,   C*q_dot(:),  w,'FaceColor',th.grn, 'EdgeColor','k','DisplayName','C\dotq (Coriolis)');
    b3 = bar(ax,x+w, G,           w,'FaceColor',th.red, 'EdgeColor','k','DisplayName','G (gravity)');
    yline(ax,0,'k-','LineWidth',1.5);
    style_axis(ax, th, 'Torque  (N.m)');
    title(ax,'Torque Component Breakdown per Joint','FontSize',13,'FontWeight','bold','Interpreter','none');
    legend(ax,[b1 b2 b3],'Location','southeast','FontSize',10,'Interpreter','none','Box','on');
end

function plot_inertia_heatmap(M, th)
    figure('Name','Inertia Matrix','Position',[60 -20 700 560],'Color','w');
    ax = axes('Position',[0.12 0.12 0.72 0.78]);
    imagesc(ax,M);
    colormap(ax,[linspace(0.95,0.12,128)', linspace(0.97,0.35,128)', linspace(1.00,0.75,128)']);
    cb = colorbar(ax); cb.Label.String = 'kg.m^2'; cb.FontSize = 10;
    set(ax,'XTick',1:6,'YTick',1:6,'XTickLabel',{'J1','J2','J3','J4','J5','J6'}, ...
           'YTickLabel',{'J1','J2','J3','J4','J5','J6'},'FontSize',12,'TickLength',[0 0]);
    axis(ax,'square');
    title(ax,'Inertia Matrix  M(q)  [kg.m^2]','FontSize',13,'FontWeight','bold','Interpreter','none');
    Mmax = max(M(:));
    for i = 1:6
        for j = 1:6
            tc = 'k'; if M(i,j)/Mmax > 0.55, tc = 'w'; end
            text(ax,j,i,sprintf('%.4f',M(i,j)),'HorizontalAlignment','center', ...
                 'VerticalAlignment','middle','FontSize',9,'FontWeight','bold','Color',tc);
        end
    end
end

function plot_robot_posture(T, q_deg, th)
    figure('Name','Robot Posture','Position',[800 -20 700 560],'Color','w');
    ax = axes; hold(ax,'on');
    pts = zeros(3,7); for i = 1:7, pts(:,i) = T{i}(1:3,4); end
    plot3(ax,pts(1,:),pts(2,:),pts(3,:),'-','Color',th.ink,'LineWidth',4.5);
    plot3(ax,pts(1,:),pts(2,:),pts(3,:),'o','MarkerSize',10, ...
          'MarkerFaceColor',th.red,'MarkerEdgeColor','k','LineWidth',1.2);
    plot3(ax,pts(1,7),pts(2,7),pts(3,7),'p','MarkerSize',16, ...
          'MarkerFaceColor',th.gold,'MarkerEdgeColor','k');
    sc = 0.05;
    draw_frame(ax, eye(3), [0;0;0], sc);            % base frame
    draw_frame(ax, T{7}(1:3,1:3), T{7}(1:3,4), sc); % end-effector frame
    grid(ax,'on'); axis(ax,'equal'); box(ax,'on');
    xlabel(ax,'X (m)','FontWeight','bold'); ylabel(ax,'Y (m)','FontWeight','bold'); zlabel(ax,'Z (m)','FontWeight','bold');
    title(ax,sprintf('Robot Posture   [%.1f  %.1f  %.1f  %.1f  %.1f  %.1f] deg',q_deg), ...
          'FontSize',12,'FontWeight','bold','Interpreter','none');
    view(ax,135,25); set(ax,'FontSize',11,'GridAlpha',0.25);
end

function draw_frame(ax, R, o, sc)
    cols = {'r','g','b'};
    for k = 1:3
        quiver3(ax,o(1),o(2),o(3),R(1,k)*sc,R(2,k)*sc,R(3,k)*sc, ...
                cols{k},'LineWidth',2,'MaxHeadSize',0.8,'AutoScale','off');
    end
end

function style_axis(ax, th, ylab)
    set(ax,'XTick',1:6,'XTickLabel',{'Joint 1','Joint 2','Joint 3','Joint 4','Joint 5','Joint 6'}, ...
           'FontSize',th.font,'Box','on','GridAlpha',0.3);
    xlim(ax,[0.4 6.6]); grid(ax,'on');
    ylabel(ax,ylab,'FontSize',12,'FontWeight','bold');
end
