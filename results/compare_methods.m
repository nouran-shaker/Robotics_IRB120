%% =====================================================================
%  COMPARE  -  your analytical modules  vs  MATLAB toolbox
%  ---------------------------------------------------------------------
%  Builds a rigidBodyTree that matches your jacobian.m / irb120_dynamics.m
%  EXACTLY (same DH, same masses/COMs/inertias, metres), then compares:
%     * geometric Jacobian   (jacobian.m        vs geometricJacobian)
%     * inertia matrix M(q)  (irb120_dynamics   vs massMatrix)
%     * gravity torque G(q)  (irb120_dynamics   vs gravityTorque)
%     * Coriolis  C*qd       (irb120_dynamics   vs velocityProduct)
%     * full tau             (irb120_dynamics   vs inverseDynamics)
%
%  THREE alignment points handled here:
%   1) Jacobian row order: yours = [linear; angular], toolbox = [angular;
%      linear]  ->  toolbox J is row-swapped before comparing.
%   2) Units: your scripts are in METRES, so the tree is built in metres.
%   3) Inertia: your p.I is about the COM; rigidBody.Inertia is about the
%      body-frame ORIGIN  ->  parallel-axis shift applied per link.
%
%  Requires jacobian.m and irb120_dynamics.m on the path.
%  (jacobian.m runs clc internally, so results are collected and printed
%   once at the very end.)
% =====================================================================

clc; clear; close all;

P = get_params_local();          % same numbers as irb120_dynamics get_params
robot  = build_tree(P);
eeName = 'link6';

% test configurations (deg) - all within joint limits
Qdeg = [ 10 -30  40  20  30 -15;
        -45  30 -20  60  40 -30;
         20 -60  30 -40  50  25;
         30  20  10  15 -25  35 ];
nQ = size(Qdeg,1);

% storage (printed after the loop, because jacobian.m calls clc)
jac = struct('dJ',[],'mA',[],'mT',[]);
dyn = struct('dM',[],'dG',[],'dC',[],'dTau',[]);

jOpts = struct('Plots',false,'Verbose',false);

for c = 1:nQ
    q   = deg2rad(Qdeg(c,:));
    qd  = 0.5*[ 1 -1 0.5 -0.5 1 -1];
    qdd = 0.3*[-1  1 -0.5 0.5 -1 1];

    % ---------- JACOBIAN ----------
    Tf = fk_frames(q, P);                          % {1..7} metres
    Ja = jacobian(q, Tf, qd(:), jOpts);            % yours: [linear; angular]
    Jt = geometricJacobian(robot, q, eeName);      % toolbox: [angular; linear]
    Jt_la = [Jt(4:6,:); Jt(1:3,:)];                % -> [linear; angular]
    jac.dJ(c) = max(abs(Ja(:) - Jt_la(:)));
    jac.mA(c) = sqrt(max(det(Ja*Ja'),0));
    jac.mT(c) = sqrt(max(det(Jt*Jt'),0));

    % ---------- DYNAMICS ----------
    r = irb120_dynamics(q, qd, qdd, 'Plots',false,'Verbose',false,'Export',false);
    Mt   = massMatrix(robot, q);
    Gt   = gravityTorque(robot, q).';
    Cqdt = velocityProduct(robot, q, qd).';
    tauT = inverseDynamics(robot, q, qd, qdd).';

    dyn.dM(c)   = max(abs(r.M(:)  - Mt(:)));
    dyn.dG(c)   = max(abs(r.G     - Gt));
    dyn.dC(c)   = max(abs(r.C*qd(:) - Cqdt));
    dyn.dTau(c) = max(abs(r.tau   - tauT));
end

%% ---------------- print everything (after clc-ing calls) -------------
fprintf('\n======================================================\n');
fprintf('   JACOBIAN   yours vs geometricJacobian (row-aligned)\n');
fprintf('======================================================\n');
fprintf('  %-7s | maxAbsDiff |  manip(yours)   manip(toolbox)\n','config');
fprintf('  %s\n', repmat('-',1,55));
for c = 1:nQ
    fprintf('  q%-6d | %10.3e | %12.5f   %12.5f\n', c, jac.dJ(c), jac.mA(c), jac.mT(c));
end

fprintf('\n======================================================\n');
fprintf('   DYNAMICS   yours vs toolbox   (max abs diff)\n');
fprintf('======================================================\n');
fprintf('  %-7s |    M(q)    |    G(q)    |   C*qd    |    tau\n','config');
fprintf('  %s\n', repmat('-',1,58));
for c = 1:nQ
    fprintf('  q%-6d | %9.2e | %9.2e | %9.2e | %9.2e\n', ...
            c, dyn.dM(c), dyn.dG(c), dyn.dC(c), dyn.dTau(c));
end

fprintf('\n======================================================\n');
fprintf('  Expected: J, M, G match ~1e-9 (closed form both sides).\n');
fprintf('  C*qd and tau match ~1e-6 (your C is central-difference,\n');
fprintf('  toolbox C is exact). Larger M/G diff => check the inertia\n');
fprintf('  parallel-axis shift or COM/mass values.\n');
fprintf('======================================================\n');


%% =====================================================================
%  LOCAL FUNCTIONS
% =====================================================================
function P = get_params_local()
    % --- MUST mirror irb120_dynamics get_params() ---
    P.d1=0.290; P.a2=0.270; P.a3=0.070; P.d4=0.302; P.d6=0.072;   % m
    P.g  = [0 0 -9.81];
    P.m  = [7.0 6.5 4.0 3.5 2.5 1.5];
    P.com = [0.000 0.000 0.145;
             0.135 0.000 0.000;
             0.035 0.000 0.000;
             0.000 0.000 0.151;
             0.000 0.000 0.000;
             0.000 0.000 0.036];
    rad = [0.090 0.050 0.040 0.030 0.025 0.040];
    len = [P.d1 P.a2 P.a3 P.d4 0 P.d6];
    P.I = cell(1,6);
    for i = 1:6
        m=P.m(i); rr=rad(i); L=len(i);
        if i==5
            Ixx=(2/5)*m*rr^2; Iyy=Ixx; Izz=Ixx;
        else
            Ixx=(1/12)*m*(3*rr^2+L^2); Iyy=Ixx; Izz=(1/2)*m*rr^2;
        end
        P.I{i} = diag([Ixx Iyy Izz]);
    end
end

function robot = build_tree(P)
    dhRows = [    0   -pi/2   P.d1   0;
               P.a2     0      0     0;
               P.a3  -pi/2     0     0;
                  0   -pi/2   P.d4   0;
                  0    pi/2     0     0;
                  0      0    P.d6   0];
    robot = rigidBodyTree('DataFormat','row');
    robot.Gravity = P.g;
    parent = robot.BaseName;
    for i = 1:6
        b = rigidBody(sprintf('link%d',i));
        j = rigidBodyJoint(sprintf('joint%d',i),'revolute');
        setFixedTransform(j, dhRows(i,:), 'dh');
        d  = P.com(i,:).';
        Ic = P.I{i};
        Io = Ic + P.m(i)*((d.'*d)*eye(3) - d*d.');   % COM -> body origin
        b.Mass = P.m(i);
        b.CenterOfMass = P.com(i,:);
        b.Inertia = [Io(1,1) Io(2,2) Io(3,3) Io(2,3) Io(1,3) Io(1,2)];
        b.Joint = j;
        addBody(robot, b, parent);
        parent = b.Name;
    end
end

function T = fk_frames(q, P)
    T = cell(7,1); T{1}=eye(4);
    T{2}=T{1}*dh_mat(q(1),-pi/2,0,P.d1);
    T{3}=T{2}*dh_mat(q(2),0,P.a2,0);
    T{4}=T{3}*dh_mat(q(3),-pi/2,P.a3,0);
    T{5}=T{4}*dh_mat(q(4),-pi/2,0,P.d4);
    T{6}=T{5}*dh_mat(q(5),pi/2,0,0);
    T{7}=T{6}*dh_mat(q(6),0,0,P.d6);
end

function A = dh_mat(theta,alpha,a,d)
    ct=cos(theta);st=sin(theta);ca=cos(alpha);sa=sin(alpha);
    A=[ct,-ca*st,sa*st,a*ct; st,ca*ct,-sa*ct,a*st; 0,sa,ca,d; 0,0,0,1];
end