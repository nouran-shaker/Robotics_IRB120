clc; clear workspace;
DOF2_Arm.DataFormat = "column"; % important
N = 15000; % number of samples
workspace = zeros(N,3);
eeName = DOF2_Arm.BodyNames{end}; % usually last body is end-effector
for i = 1:N
config = randomConfiguration(DOF2_Arm); % respects joint limits
T = getTransform(DOF2_Arm,config,eeName);
workspace(i,:) = T(1:3,4)';
end
figure;
scatter3(workspace(:,1),workspace(:,2),workspace(:,3),5,'filled');
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Workspace Point Cloud');
view(2); % for planar RR robot

figure
scatter3(workspace(:,1),workspace(:,2),workspace(:,3),5,'filled')
axis equal
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Robot Workspace (3D)')
view(3)