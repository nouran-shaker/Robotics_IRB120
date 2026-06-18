function workspace_analysis(opts)
%WORKSPACE_ANALYSIS  Reachable workspace of the ABB IRB120 by Monte-Carlo
%   forward-kinematics sampling, with a voxel volume estimate.
%
%   workspace_analysis()       defaults (~400k sample points)
%   workspace_analysis(opts)   opts.N / opts.Plots / opts.Voxel / opts.DexThresh
%
%   NOTE on "configurations": the reachable workspace is configuration-
%   INDEPENDENT. The 4/8 elbow-wrist-shoulder branches are inverse-kinematics
%   solutions for reaching ONE given pose; the workspace is the set of ALL
%   reachable tip positions, so it carries no solution split. The cloud is
%   therefore shown as a single set, coloured by height.

if nargin < 1, opts = struct(); end
if ~isfield(opts,'N'),         opts.N         = 400000; end
if ~isfield(opts,'Plots'),     opts.Plots     = true;   end
if ~isfield(opts,'Voxel'),     opts.Voxel     = 0.040;  end   % 40 mm
if ~isfield(opts,'DexThresh'), opts.DexThresh = 0.50;   end

fprintf('====================================================\n');
fprintf('   ABB IRB120 Workspace Analysis  (reachable)\n');
fprintf('====================================================\n\n');

%% ---- robot parameters (m) ----
params.d1=0.290; params.a2=0.270; params.a3=0.070; params.d4=0.302; params.d6=0.072;
limits_deg = [-165 165; -110 110; -110 70; -160 160; -120 120; -400 400];
limits_rad = deg2rad(limits_deg);

%% ---- sample joints, record tip position AND approach direction ----
N = opts.N;
fprintf('Sampling %d points (joints 1-5; joint 6 leaves the tip pose unchanged)...\n', N);
q = zeros(N,6);
for j = 1:5
    q(:,j) = limits_rad(j,1) + (limits_rad(j,2)-limits_rad(j,1)).*rand(N,1);
end
allP = zeros(N,3); allA = zeros(N,3);
for i = 1:N
    T = fk_position(q(i,:), params);
    allP(i,:) = T(1:3,4)';
    allA(i,:) = T(1:3,3)';          % approach vector a = tool z-axis
end
fprintf('Point cloud complete.\n');

%% ---- reachable workspace volume (voxel occupancy) ----
vsize = opts.Voxel;
key   = floor(allP / vsize);
[ukey,~,ic] = unique(key,'rows');
nVox      = size(ukey,1);
reach_vol = nVox * vsize^3;

%% ---- dexterous workspace VOLUME (orientation coverage per voxel) ----
% kept as a numeric metric (Part B deliverable); no separate figure.
[dx,dy,dz] = ndgrid(-1:1,-1:1,-1:1);
dirs = [dx(:) dy(:) dz(:)];
dirs(all(dirs==0,2),:) = [];
dirs = dirs ./ sqrt(sum(dirs.^2,2));            % 26 direction bins
nDir = size(dirs,1);
[~,binIdx] = max(allA * dirs', [], 2);
occ        = sparse(ic, binIdx, 1, nVox, nDir);
cov_frac   = full(sum(occ > 0, 2)) / nDir;
dex_mask   = cov_frac >= opts.DexThresh;
dex_vol    = sum(dex_mask) * vsize^3;
vox_ctr    = (ukey + 0.5) * vsize;

%% ============================================================
%  FIGURES   (1) 3D reachable cloud   (2) top + front projections
%% ============================================================
if opts.Plots
    % Fig 1: single reachable cloud, coloured by height
    figure('Name','IRB120 Reachable Workspace','Position',[50 50 900 750],'Color','w');
    scatter3(allP(:,1),allP(:,2),allP(:,3),1,allP(:,3),'.');
    hold on; plot3(0,0,0,'k^','MarkerSize',12,'MarkerFaceColor','k');
    grid on; axis equal; box on; view(135,25); colormap(jet);
    cb = colorbar; cb.Label.String = 'Z (m)';
    xlabel('X (m)','FontWeight','bold'); ylabel('Y (m)','FontWeight','bold'); zlabel('Z (m)','FontWeight','bold');
    title({'ABB IRB120 - Reachable Workspace', sprintf('%d sample points', N)}, ...
          'FontSize',14,'FontWeight','bold');

    % Fig 2: height-coloured orthographic projections (top + front)
    figure('Name','IRB120 Workspace - Projections','Position',[200 80 1150 520],'Color','w');
    subplot(1,2,1);
    scatter(allP(:,1),allP(:,2),1,allP(:,3),'.');
    axis equal; box on; grid on; colormap(jet);
    cb = colorbar; cb.Label.String = 'Z (m)';
    xlabel('X (m)','FontWeight','bold'); ylabel('Y (m)','FontWeight','bold');
    title('Top view  (perpendicular to Z)','FontSize',13,'FontWeight','bold');
    subplot(1,2,2);
    scatter(allP(:,1),allP(:,3),1,allP(:,3),'.');
    axis equal; box on; grid on; colormap(jet);
    cb = colorbar; cb.Label.String = 'Z (m)';
    xlabel('X (m)','FontWeight','bold'); ylabel('Z (m)','FontWeight','bold');
    title('Front view  (perpendicular to Y)','FontSize',13,'FontWeight','bold');
end

%% ============================================================
%  STATISTICS
%% ============================================================
fprintf('\n====================================================\n');
fprintf('   WORKSPACE STATISTICS\n');
fprintf('====================================================\n');
xr=[min(allP(:,1)) max(allP(:,1))]; yr=[min(allP(:,2)) max(allP(:,2))]; zr=[min(allP(:,3)) max(allP(:,3))];
fprintf('  X: [%.3f, %.3f] m   span %.3f m\n', xr(1),xr(2),diff(xr));
fprintf('  Y: [%.3f, %.3f] m   span %.3f m\n', yr(1),yr(2),diff(yr));
fprintf('  Z: [%.3f, %.3f] m   span %.3f m\n', zr(1),zr(2),diff(zr));
r_max    = max(sqrt(allP(:,1).^2 + allP(:,2).^2 + (allP(:,3)-params.d1).^2));
r_theory = params.a2 + sqrt(params.a3^2+params.d4^2) + params.d6;
fprintf('  Max reach (sampled):  %.3f m (%.1f mm)\n', r_max, r_max*1000);
fprintf('  Max reach (theory):   %.3f m (%.1f mm)\n', r_theory, r_theory*1000);
fprintf('  --- volumes (voxel %.0f mm) ---\n', vsize*1000);
fprintf('  Reachable workspace:  %.4f m^3  (%d voxels)\n', reach_vol, nVox);
fprintf('  Dexterous workspace:  %.4f m^3  (%d voxels, >= %.0f%% orient. coverage)\n', ...
        dex_vol, sum(dex_mask), 100*opts.DexThresh);
if reach_vol > 0
    fprintf('  Dexterous / reachable: %.1f %%\n', 100*dex_vol/reach_vol);
end
fprintf('====================================================\n');

%% ---- export ----
assignin('base','workspace_pts',allP);
assignin('base','workspace_app',allA);
assignin('base','workspace_reach_vol',reach_vol);
assignin('base','workspace_dex_vol',dex_vol);
assignin('base','workspace_dex_centers',vox_ctr(dex_mask,:));
fprintf('\nExported: workspace_pts, workspace_app, workspace_reach_vol,\n');
fprintf('          workspace_dex_vol, workspace_dex_centers\n');
fprintf('====================================================\n');
end


%% ============================================================
%  LOCAL FK
%% ============================================================
function T = fk_position(q, DH)
    T = dh_mat(q(1),-pi/2,0,DH.d1) * dh_mat(q(2),0,DH.a2,0) * ...
        dh_mat(q(3),-pi/2,DH.a3,0) * dh_mat(q(4),-pi/2,0,DH.d4) * ...
        dh_mat(q(5),pi/2,0,0) * dh_mat(q(6),0,0,DH.d6);
end

function A = dh_mat(theta,alpha,a,d)
    ct=cos(theta); st=sin(theta); ca=cos(alpha); sa=sin(alpha);
    A = [ct,-ca*st,sa*st,a*ct; st,ca*ct,-sa*ct,a*st; 0,sa,ca,d; 0,0,0,1];
end
