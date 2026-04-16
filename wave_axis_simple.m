%% wave_axis_simple.m
% Minimal single-file version of the wave-axis analysis.
% Loads a boundary_flows results .mat, asks you to click the nucleus, and
% tells you whether the wave axis points at it. Expect ~60 lines.
%
% Set `resultsFile` in the workspace first, e.g.
%   resultsFile = fullfile(outDir, 'tangential_linear_interp_results.mat');
% then just run this script.

%% 1. Load
assert(exist('resultsFile','var')==1 && exist(resultsFile,'file')==2, ...
    'Set resultsFile to a .mat produced by one of the boundary_flows scripts.');
S = load(resultsFile);
nF = size(S.Vtheta_kymo, 1);

%% 2. Click the nucleus on the first frame
if isfield(S, 'base_dir')
    dlist = dir(fullfile(S.base_dir, '*State1*'));
    dlist = dlist(arrayfun(@(d) ~d.isdir && d.name(1) ~= '.', dlist));
    I0 = double(imread(fullfile(dlist(1).folder, dlist(1).name)));
else
    I0 = zeros(1024);  % blank fallback
end
figure; imshow(I0, []); hold on;
cxy = S.centroidXY(find(S.qcFlag,1), :);
plot(cxy(1), cxy(2), 'c+', 'MarkerSize', 14, 'LineWidth', 2);
title('Click the nucleus center');
[xn, yn] = ginput(1); close;

%% 3. Per-frame wave axis from velocity antisymmetry + nucleus direction
phi_vel   = nan(nF, 1);
theta_nuc = nan(nF, 1);
haveRaw   = isfield(S, 'Vtheta_raw') && isfield(S, 'Theta_raw');

for fr = 1:nF
    if ~S.qcFlag(fr), continue; end
    if haveRaw
        v = S.Vtheta_raw{fr}(:);   th = S.Theta_raw{fr}(:);
    else
        v = S.Vtheta_kymo(fr,:).'; th = S.thetaCenters(1:numel(v)).';
    end
    ok = isfinite(v) & isfinite(th);
    if nnz(ok) < 4, continue; end
    phi_vel(fr) = mod(pi/2 - angle(sum(v(ok) .* exp(-1i*th(ok)))), 2*pi);

    c = S.centroidXY(fr,:);
    theta_nuc(fr) = mod(atan2(yn - c(2), xn - c(1)), 2*pi);
end

dPhi = angle(exp(1i*(phi_vel - theta_nuc)));   % wrapped to [-pi, pi]

%% 4. Circular statistics over all QC frames
d = dPhi(isfinite(dPhi));
r_mean = mean(exp(1i*d));
circ_mu   = angle(r_mean);
circ_R    = abs(r_mean);
ray_p     = exp(-numel(d) * circ_R^2);
fprintf('\nWave axis vs nucleus direction (%d frames):\n', numel(d));
fprintf('  circular mean Dphi = %.1f deg\n', rad2deg(circ_mu));
fprintf('  resultant length R = %.3f  (0 = random, 1 = perfectly aligned)\n', circ_R);
fprintf('  Rayleigh p ~ %.3g  (<0.05 = significant preferred direction)\n', ray_p);

%% 5. One summary figure
figure('Position', [100 100 1200 500]);
subplot(1,2,1);
imagesc(rad2deg(S.thetaCenters), S.time_min, S.Vtheta_kymo); hold on;
plot(rad2deg(phi_vel),   S.time_min, 'w-',  'LineWidth', 1.5);
plot(rad2deg(theta_nuc), S.time_min, 'm--', 'LineWidth', 1.5);
xlabel('\theta (deg)'); ylabel('Time (min)');
title('v_\theta(\theta,t): \phi_{vel} (white), \theta_{nuc} (magenta)');
colorbar;

subplot(1,2,2);
polarhistogram(d, 24, 'FaceColor', [.2 .5 .8]);
title(sprintf('\\Delta\\phi distribution\nmean=%.1f\\circ, R=%.2f, p=%.2g', ...
    rad2deg(circ_mu), circ_R, ray_p));
