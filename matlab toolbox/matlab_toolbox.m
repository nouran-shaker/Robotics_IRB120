%% =====================================================================
%  ABB IRB120 - TOOLBOX METHOD, built on the ANALYTICAL DH FRAMES
%  The matched DH 'robot' does ALL analysis (FK/IK/Jacobian/dynamics/
%  comparison). The animation uses the loadrobot() MESH model 'robotViz',
%  driven by its OWN IK on the SAME Cartesian trajectory - because the two
%  models use different frame conventions, so a given joint-angle vector is
%  NOT the same physical pose in both. Solving robotViz's own IK makes the
%  mesh tool tip follow the trajectory line correctly.
%
%  Requires inverse_kinematics.m on the path. Comment out its top line
%  "clc; close all;" so the comparison loop doesn't wipe the console.
% =====================================================================

clc; clear; close all;

%% ==============================
% 1. BUILD ROBOT FROM ANALYTICAL DH  (replaces loadrobot)
% ===============================
DH.d1 = 290; DH.a2 = 270; DH.a3 = 70; DH.d4 = 302; DH.d6 = 72;   % mm

% setFixedTransform 'dh' wants rows of [a  alpha  d  theta]
dhRows = [    0   -pi/2   DH.d1   0;
           DH.a2     0      0     0;
           DH.a3  -pi/2     0     0;
              0   -pi/2   DH.d4   0;
              0    pi/2     0     0;
              0      0    DH.d6   0];
limDeg = [-165 165; -110 110; -110 70; -160 160; -120 120; -400 400];

robot = rigidBodyTree('DataFormat','row');
robot.Gravity = [0 0 -9810];            % mm/s^2
eeName = 'link6';                        % == analytical tool tip
parent = robot.BaseName;
for i = 1:6
    b = rigidBody(sprintf('link%d',i));
    j = rigidBodyJoint(sprintf('joint%d',i),'revolute');
    setFixedTransform(j, dhRows(i,:), 'dh');
    j.PositionLimits = deg2rad(limDeg(i,:));
    b.Mass = 3.0; b.CenterOfMass = [0 0 0];        % placeholder inertia so
    b.Inertia = [1e4 1e4 1e4 0 0 0];               % inverseDynamics runs
    b.Joint = j;
    addBody(robot, b, parent);
    parent = b.Name;
end

% Visual twin with meshes (metres) — used ONLY for the animation display
robotViz = loadrobot('abbIrb120','DataFormat','row');

figure(1);
show(robotViz,'Visuals','on');
title('ABB IRB120');

%% ==============================
% 2. HOME CONFIGURATION
% ===============================
q_home = homeConfiguration(robot);

%% ==============================
% 3. TOOL ORIENTATION (VERTICAL)
% ===============================
R_down = eul2tform([pi 0 0]); % tool pointing down

%% ==============================
% 4. IK SETUP
% ===============================
ik = inverseKinematics('RigidBodyTree',robot);
weights = [0.5 0.5 0.5 1 1 1];

%% ==============================
% 5. PICK & PLACE POSITIONS  (mm)
% ===============================
pickPos = [300  200  50];
pickA   = [300  200 150];

placePos = [200 -300  50];
placeA   = [200 -300 150];

%% ==============================
% 6. CARTESIAN TRAJECTORY (FIXED)
% ===============================
n = 40;

T1 = transformtraj(trvec2tform(pickA)*R_down,...
                   trvec2tform(pickPos)*R_down,[0 1],linspace(0,1,n));

T2 = transformtraj(trvec2tform(pickPos)*R_down,...
                   trvec2tform(pickA)*R_down,[0 1],linspace(0,1,n));

T3 = transformtraj(trvec2tform(pickA)*R_down,...
                   trvec2tform(placeA)*R_down,[0 1],linspace(0,1,n));

T4 = transformtraj(trvec2tform(placeA)*R_down,...
                   trvec2tform(placePos)*R_down,[0 1],linspace(0,1,n));

T5 = transformtraj(trvec2tform(placePos)*R_down,...
                   trvec2tform(placeA)*R_down,[0 1],linspace(0,1,n));

T_all = cat(3,T1,T2,T3,T4,T5);

%% ==============================
% 7. IK TRAJECTORY  (matched DH robot - for analysis)
% ===============================
qTraj = zeros(6,size(T_all,3));

for i = 1:size(T_all,3)
    [qTraj(:,i),~] = ik(eeName,T_all(:,:,i),weights,q_home);
end

%% ==============================
% 8. END EFFECTOR PATH  (matched DH robot, mm)
% ===============================
eePath = zeros(size(qTraj,2),3);

for i = 1:size(qTraj,2)
    T = getTransform(robot,qTraj(:,i)',eeName);
    eePath(i,:) = tform2trvec(T);
end

figure(2);
plot3(eePath(:,1),eePath(:,2),eePath(:,3),'b','LineWidth',2);
grid on; title('End Effector Path');

%% ==============================
% 9. WORKSPACE
% ===============================
N = 500;
ws = zeros(N,3);

for i = 1:N
    q_rand = randomConfiguration(robot);
    T = getTransform(robot,q_rand,eeName);
    ws(i,:) = tform2trvec(T);
end

figure(3);
scatter3(ws(:,1),ws(:,2),ws(:,3),'.');
title('Workspace'); grid on;

%% ==============================
% 10. JACOBIAN + SINGULARITY
% ===============================
q_test = [0.2 -0.5 0.3 0.1 0.2 -0.1];
J = geometricJacobian(robot,q_test,eeName);
detJ = det(J*J');

%% ==============================
% 11. INVERSE DYNAMICS
% ===============================
qd = gradient(qTraj')';
qdd = gradient(qd')';

tau = zeros(size(qTraj'));

for i = 1:length(qTraj)
    tau(i,:) = inverseDynamics(robot,...
        qTraj(:,i)',qd(:,i)',qdd(:,i)');
end

figure(4);
plot(tau);
title('Joint Torques'); legend('J1','J2','J3','J4','J5','J6');

%% ==============================
% 12. TRAJECTORY ANALYSIS
% ===============================
figure(5);
subplot(3,1,1); plot(qTraj'); title('Position');
subplot(3,1,2); plot(qd'); title('Velocity');
subplot(3,1,3); plot(qdd'); title('Acceleration');

%% ==============================
% 13. JOINT LIMITS
% ===============================
disp('Joint Limits:');
for i = 1:length(robot.Bodies)
    joint = robot.Bodies{i}.Joint;
    if strcmp(joint.Type,'revolute')
        disp(joint.PositionLimits);
    end
end

%% ==============================
% 14. DH (TRANSFORMS)
% ===============================
disp('Joint Transforms:');
for i = 1:length(robot.Bodies)
    disp(robot.Bodies{i}.Joint.JointToParentTransform);
end

%% ==============================
% 14b. METHOD COMPARISON  (closed-form vs toolbox, same DH model)
%      Placed here because the animation loop below never returns.
% ===============================
fprintf('\n=====================================================\n');
fprintf('   CLOSED-FORM  vs  TOOLBOX   (identical DH model)\n');
fprintf('=====================================================\n');

% FK sanity: toolbox vs analytical FK
qchk = deg2rad([20 -30 25 40 -50 15]);
T_tb = getTransform(robot,qchk,eeName);
T_an = fk_analytic(qchk,DH);
fprintf('FK match  |  max |toolbox - analytic| = %.3e mm\n', max(abs(T_tb(:)-T_an(:))));

testPoses = { ...
 [ 0.1283 0.1969 0.9720 455.98; 0.8884 0.4128 -0.2009 368.52; -0.4408 0.8893 -0.1219 394.35; 0 0 0 1], ...
 [ 0.0682 0.6780 -0.7319 299.04; 0.5376 0.5930 0.5994 445.68;  0.8404 -0.4343 -0.3240 41.87; 0 0 0 1] };

opts.Plot = false; opts.Interactive = false; opts.CloudPts = 0;   % skip cloud -> fast
for tp = 1:numel(testPoses)
    T = testPoses{tp};
    fprintf('\n##### Test pose %d #####\n', tp);
    cf = inverse_kinematics(T, opts);          % analytical, all 8 branches
    if ~isfield(cf,'is_valid') || ~any(cf.is_valid)
        fprintf('  No valid analytical solution.\n'); continue;
    end
    vIdx = find(cf.is_valid);

    fprintf('\n  [A] seed toolbox with each analytical branch (expect ~0):\n');
    fprintf('  %-26s | maxAngDiff(deg) | posErr(mm) | rotErr\n','branch');
    fprintf('  %s\n', repmat('-',1,70));
    for m = 1:numel(vIdx)
        k=vIdx(m); q_an=cf.theta_rad(k,:);
        q_tb = ik(eeName,T,weights,q_an);
        dAng = max(abs(rad2deg(wrapPi(q_tb-q_an))));
        Ttb  = getTransform(robot,q_tb,eeName);
        pErr = norm(Ttb(1:3,4)-T(1:3,4));
        rErr = norm(Ttb(1:3,1:3)-T(1:3,1:3),'fro');
        fprintf('  %-26s | %14.4f | %10.3e | %.2e\n', cf.config_names{k}, dAng, pErr, rErr);
    end

    q_ind = ik(eeName,T,weights,q_home);
    d = inf(1,numel(vIdx));
    for m=1:numel(vIdx), d(m)=max(abs(rad2deg(wrapPi(q_ind-cf.theta_rad(vIdx(m),:))))); end
    [dmin,mb]=min(d); kb=vIdx(mb);
    fprintf('\n  [B] toolbox solved from home -> branch "%s" (diff %.4f deg)\n', cf.config_names{kb}, dmin);
    fprintf('      toolbox q (deg): [%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f]\n', rad2deg(q_ind));
    fprintf('      analyt. q (deg): [%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f]\n', cf.theta_deg(kb,:));
end
fprintf('\n  Analytical returns ALL 8 closed-form solutions; toolbox returns\n');
fprintf('  ONE (depends on the seed). [A] ~0 confirms identical frames.\n');
fprintf('=====================================================\n');

%% ==============================
% RESULTS WINDOW (MAIN DISPLAY)
% ===============================
figResults = figure('Name','Robot Analysis Results',...
    'NumberTitle','off',...
    'Position',[50 100 750 500]);

T_fk = getTransform(robot,q_home,eeName);
[qSol,~] = ik(eeName,trvec2tform([300 200 200])*R_down,weights,q_home);
J = geometricJacobian(robot,q_home,eeName);
detJ = det(J*J');

% Titles
uicontrol('Parent',figResults,'Style','text',...
    'Position',[50 460 200 20],...
    'String','Forward Kinematics','FontWeight','bold');

uicontrol('Parent',figResults,'Style','text',...
    'Position',[350 460 200 20],...
    'String','Inverse Kinematics','FontWeight','bold');

uicontrol('Parent',figResults,'Style','text',...
    'Position',[50 240 200 20],...
    'String','Jacobian','FontWeight','bold');

uicontrol('Parent',figResults,'Style','text',...
    'Position',[600 240 150 20],...
    'String','Det(JJ^T)','FontWeight','bold');

% Tables
uitable('Parent',figResults,'Data',T_fk,...
    'Position',[50 280 250 150]);

uitable('Parent',figResults,'Data',qSol,...
    'Position',[350 320 250 60]);

uitable('Parent',figResults,'Data',J,...
    'Position',[50 20 500 200]);

uitable('Parent',figResults,'Data',detJ,...
    'Position',[600 120 120 60]);

drawnow;
uiwait(msgbox('Review results, then press OK to start simulation','Ready'));

%% ==============================
% ANIMATION PREP  (solve the VISUAL model's own IK on the same path)
%   robotViz uses a different frame convention + metres, so it must solve
%   its own IK; the mesh then tracks the trajectory line exactly.
% ===============================
ikViz    = inverseKinematics('RigidBodyTree',robotViz);
qHomeViz = homeConfiguration(robotViz);

T_all_m = T_all;
T_all_m(1:3,4,:) = T_all(1:3,4,:)/1000;        % mm -> m on translation only

qTrajViz = zeros(6,size(T_all_m,3));
seed = qHomeViz;
for i = 1:size(T_all_m,3)
    [qTrajViz(:,i),~] = ikViz('tool0',T_all_m(:,:,i),weights,seed);
    seed = qTrajViz(:,i)';                      % warm-start
end

eePathViz = zeros(size(qTrajViz,2),3);          % loadrobot tool path (m)
for i = 1:size(qTrajViz,2)
    Tv = getTransform(robotViz,qTrajViz(:,i)','tool0');
    eePathViz(i,:) = tform2trvec(Tv);
end

%% ==============================
% ANIMATION WINDOW  (solid mesh model on the trajectory)
% ===============================
figAnim = figure('Name','Robot Simulation',...
    'NumberTitle','off',...
    'Position',[850 100 700 500]);

axes;
hold on; grid on;
view(135,25);
axis equal;

plot3(eePathViz(:,1),eePathViz(:,2),eePathViz(:,3),'r--','LineWidth',1.5);

trail = animatedline('Color','b','LineWidth',2);

set(figAnim,'Renderer','opengl')

while true

    clearpoints(trail);

    for i = 1:size(qTrajViz,2)

        show(robotViz,qTrajViz(:,i)',...
            'PreservePlot',false,...
            'Visuals','on',...
            'Collisions','off',...
            'Frames','off');

        camlight headlight
        lighting gouraud

        addpoints(trail,eePathViz(i,1),eePathViz(i,2),eePathViz(i,3));

        drawnow limitrate;
        pause(0.02);
    end
end

%% ==============================
% LOCAL FUNCTIONS (analytical FK mirror)
% ===============================
function T = fk_analytic(q, DH)
    T = dhmat(q(1),-pi/2,0,DH.d1)*dhmat(q(2),0,DH.a2,0)* ...
        dhmat(q(3),-pi/2,DH.a3,0)*dhmat(q(4),-pi/2,0,DH.d4)* ...
        dhmat(q(5),pi/2,0,0)*dhmat(q(6),0,0,DH.d6);
end

function A = dhmat(theta,alpha,a,d)
    ct=cos(theta);st=sin(theta);ca=cos(alpha);sa=sin(alpha);
    A=[ct,-ca*st,sa*st,a*ct; st,ca*ct,-sa*ct,a*st; 0,sa,ca,d; 0,0,0,1];
end

function y = wrapPi(x)
    y = atan2(sin(x),cos(x));
end