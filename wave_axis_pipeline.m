%% wave_axis_pipeline.m
% Post-processing driver that reads tangential_kymo_results.mat (produced by
% boundary_flows.m), computes per-frame wave propagation axes from both the
% boundary shape and the tangential velocity field, prompts the user to
% click the nucleus, compares the axes to the centroid->nucleus direction,
% detects discrete wave events, and saves the enriched results back to the
% same .mat file.
%
% USAGE:
%   1. Edit resultsDir below (or set resultsDir in the base workspace).
%   2. Run the script.
%   3. A figure will appear asking you to click the nucleus on the first
%      QC-passed frame. Click once. (For moving nuclei, edit opts.clickFrames.)
%   4. Plots are saved under <resultsDir>/wave_axis/.
%
% ADDS to tangential_kymo_results.mat:
%   nucleusXY      - [nFrames x 2] nucleus (x, y) per frame (interpolated)
%   phi_shape      - [nFrames x 1] boundary-shape wave axis (rad)
%   amp_shape      - [nFrames x 1] amplitude of dominant shape mode (px)
%   kStar_shape    - [nFrames x 1] dominant Fourier mode index (1..3)
%   phi_vel        - [nFrames x 1] velocity-field wave axis (rad)
%   amp_vel        - [nFrames x 1] amplitude of dipolar velocity component
%   theta_nuc      - [nFrames x 1] atan2(yn-yc, xn-xc) direction to nucleus
%   dPhi_shape     - [nFrames x 1] wrap_pi(phi_shape - theta_nuc)
%   dPhi_vel       - [nFrames x 1] wrap_pi(phi_vel   - theta_nuc)
%   kymo_slope     - scalar: Radon-estimated wave slope on Vtheta_kymo
%                    (rad / frame; sign = direction: +=CCW, -=CW)
%   waveEvents     - struct array with fields .frame, .time_min, .phi_shape,
%                    .phi_vel, .theta_nuc, .dPhi_shape, .amp_shape
%
% ========================================================================

%% -------- CONFIG --------
% Two ways to point this script at the data:
%   1) Set `resultsFile` in the workspace to an absolute .mat path.
%   2) Set `resultsDir` and the script will auto-detect either
%      'tangential_kymo_results.mat' (from boundary_flows.m) or
%      'tangential_linear_interp_results.mat' (from
%      boundary_flows_linear_interpolant_approach.m).
if ~exist('resultsFile','var') || isempty(resultsFile)
    if ~exist('resultsDir','var') || isempty(resultsDir)
        % EDIT THIS to point to the output folder from either boundary_flows script.
        resultsDir = fullfile(pwd, 'demo_output');
    end
    candidates = { ...
        fullfile(resultsDir, 'tangential_kymo_results.mat'), ...
        fullfile(resultsDir, 'tangential_linear_interp_results.mat')};
    resultsFile = '';
    for ii = 1:numel(candidates)
        if exist(candidates{ii}, 'file') == 2
            resultsFile = candidates{ii}; break;
        end
    end
    assert(~isempty(resultsFile), ...
        ['wave_axis_pipeline: no results .mat found under %s. ' ...
         'Expected one of:\n  %s\n  %s\n' ...
         'Set `resultsFile` or `resultsDir` in the workspace.'], ...
        resultsDir, candidates{1}, candidates{2});
else
    assert(exist(resultsFile, 'file') == 2, ...
        'wave_axis_pipeline: resultsFile does not exist: %s', resultsFile);
end
resultsDir = fileparts(resultsFile);
fprintf('Using results file: %s\n', resultsFile);

outDir = fullfile(resultsDir, 'wave_axis');
if ~exist(outDir, 'dir'); mkdir(outDir); end

fprintf('=== wave_axis_pipeline ===\n');
fprintf('Loading %s ...\n', resultsFile);
S = load(resultsFile);

% --- Sanity check required fields ---
required = {'Vtheta_kymo','thetaCenters','centroidXY','time_min','qcFlag'};
for f = required
    assert(isfield(S, f{1}), ...
        'wave_axis_pipeline: field "%s" missing from %s. Re-run the boundary_flows script.', ...
        f{1}, resultsFile);
end
% visPolySeq is optional: if missing, we skip the shape-based wave axis
% and fall back to the velocity-based axis only. Re-run the updated
% boundary_flows_* script to get the full set of metrics.
haveShape = isfield(S, 'visPolySeq');
if ~haveShape
    warning(['wave_axis_pipeline: "visPolySeq" missing from results. ' ...
             'Running in VELOCITY-ONLY mode; shape-axis (phi_shape, ' ...
             'amp_shape, kStar_shape) will be NaN. Re-run the boundary_flows ' ...
             'script to get the full analysis.']);
end

Vtheta_kymo = S.Vtheta_kymo;
thetaCenters = S.thetaCenters;
centroidXY = S.centroidXY;
time_min = S.time_min;
qcFlag = S.qcFlag(:);
if haveShape
    visPolySeq = S.visPolySeq;
else
    visPolySeq = cell(size(Vtheta_kymo, 1), 1);
end
nFrames = size(Vtheta_kymo, 1);

%% -------- NUCLEUS PICKER --------
% Load existing image sequence if stored; otherwise fall back to the
% centroid-only preview (pick_nucleus handles headless fallback).
if isfield(S, 'visImageSeq')
    imgSeq = S.visImageSeq;
else
    imgSeq = cell(nFrames, 1);  % empty cells -> pick_nucleus will still show centroid
end

% Find first QC-passed frame to click on
firstQC = find(qcFlag, 1, 'first');
if isempty(firstQC); firstQC = 1; end

fprintf('Picking nucleus on frame %d (click the nucleus center)...\n', firstQC);

% Single click (static-nucleus assumption). To allow multiple clicks for a
% moving nucleus, set clickFrames below to e.g. [firstQC, nFrames/2, nFrames].
clickFrames = firstQC;

% If visImageSeq wasn't saved, try to load the first State1 image from
% base_dir (saved in results). Otherwise fall back to a blank canvas
% sized to the boundary polygon (or to the centroid).
if isempty(imgSeq{firstQC})
    loadedImg = [];
    if isfield(S, 'base_dir') && exist(S.base_dir, 'dir')
        s1_list = dir(fullfile(S.base_dir, '*State1*'));
        s1_list = s1_list(~[s1_list.isdir]);
        s1_list = s1_list(arrayfun(@(d) d.name(1) ~= '.', s1_list));
        if ~isempty(s1_list)
            idxLoad = min(firstQC, numel(s1_list));
            try
                loadedImg = double(imread(fullfile(s1_list(idxLoad).folder, ...
                                                    s1_list(idxLoad).name)));
            catch
                loadedImg = [];
            end
        end
    end
    if ~isempty(loadedImg)
        imgSeq = {loadedImg};
    elseif haveShape && ~isempty(visPolySeq{firstQC})
        warning('wave_axis_pipeline: using boundary-only preview for nucleus click.');
        poly = visPolySeq{firstQC};
        pad = 50;
        xmax = ceil(max(poly(:,1))) + pad;
        ymax = ceil(max(poly(:,2))) + pad;
        imgSeq = {zeros(ymax, xmax)};
    else
        warning('wave_axis_pipeline: using centroid-only blank canvas for nucleus click.');
        cxy = centroidXY(firstQC, :);
        side = 2 * ceil(max(cxy) + 100);
        imgSeq = {zeros(side, side)};
    end
end

opts_pick = struct('frameIdx', clickFrames, 'nFramesTotal', nFrames);
nucleusXY = pick_nucleus(imgSeq, visPolySeq, centroidXY, opts_pick);

%% -------- PER-FRAME AXES --------
fprintf('Computing per-frame wave axes (%d frames)...\n', nFrames);

phi_shape   = nan(nFrames, 1);
amp_shape   = nan(nFrames, 1);
kStar_shape = nan(nFrames, 1);
phi_vel     = nan(nFrames, 1);
amp_vel     = nan(nFrames, 1);
theta_nuc   = nan(nFrames, 1);

for fr = 1:nFrames
    if ~qcFlag(fr); continue; end
    poly = visPolySeq{fr};
    cxy  = centroidXY(fr, :);
    nxy  = nucleusXY(fr, :);

    % --- shape axis ---
    if ~isempty(poly) && all(isfinite(cxy))
        [p1, a1, k1] = wave_axis_from_shape(poly, cxy);
        phi_shape(fr)   = p1;
        amp_shape(fr)   = a1;
        kStar_shape(fr) = k1;
    end

    % --- velocity axis ---
    [p2, a2] = wave_axis_from_velocity(Vtheta_kymo(fr, :), thetaCenters);
    phi_vel(fr) = p2;
    amp_vel(fr) = a2;

    % --- nucleus direction ---
    if all(isfinite(cxy)) && all(isfinite(nxy))
        theta_nuc(fr) = mod(atan2(nxy(2) - cxy(2), nxy(1) - cxy(1)), 2*pi);
    end
end

% --- Nucleus comparison metrics (wrap to [-pi, pi]) ---
dPhi_shape = wrap_pi(phi_shape - theta_nuc);
dPhi_vel   = wrap_pi(phi_vel   - theta_nuc);

% For undirected (k*=2) axis, fold to [-pi/2, pi/2]
undir = (kStar_shape == 2);
dPhi_shape(undir) = wrap_pi_half(dPhi_shape(undir));
undir3 = (kStar_shape == 3);
dPhi_shape(undir3) = mod(dPhi_shape(undir3) + pi/3, 2*pi/3) - pi/3;

%% -------- GLOBAL KYMOGRAPH SLOPE VIA RADON --------
% Apply Radon transform to the signed kymograph and find the orientation
% of strongest diagonal structure. The result is a single signed angular
% wave velocity for the whole movie.
K = Vtheta_kymo;
K(~isfinite(K)) = 0;
% Standardize rows to emphasize traveling structure over DC bias
K = bsxfun(@minus, K, mean(K, 2, 'omitnan'));
K(~isfinite(K)) = 0;

thetas_radon = 1:179;   % degrees (Radon convention)
R = radon(K, thetas_radon);
[~, iMax] = max(var(R, 0, 1));   % strongest sinogram column
angle_deg = thetas_radon(iMax);

% Convert Radon angle to slope in (theta_rad / frame).
% Radon angle 90deg = horizontal streaks in image = no propagation.
% Angles deviating from 90 correspond to propagation with slope
% dtheta/dt = tan(90 - angle_deg) * (dTheta_per_bin / dTime_per_row).
% We return it in the natural units of the kymograph axes.
nBins   = size(Vtheta_kymo, 2);
dTheta  = 2*pi / nBins;
if numel(time_min) >= 2
    dT = (time_min(2) - time_min(1)) * 60;   % seconds per frame
else
    dT = 1;
end
slope_theta_per_sec = tand(90 - angle_deg) * (dTheta / dT);
kymo_slope = slope_theta_per_sec;

fprintf('Radon-estimated global wave slope: %.3g rad/s (%.1f deg/min)\n', ...
    kymo_slope, rad2deg(kymo_slope) * 60);

%% -------- DISCRETE WAVE EVENTS --------
% Use shape amplitude if available; otherwise fall back to velocity amplitude.
% Also track which primary metric drove detection (for downstream stats).
if haveShape && any(isfinite(amp_shape))
    amp_src = amp_shape;   primary = 'shape';
else
    amp_src = amp_vel;     primary = 'vel';
end
amp_s = amp_src;
amp_s(~isfinite(amp_s)) = 0;
win = max(3, round(5));    % 5-frame moving mean (tweak as needed)
amp_sm = movmean(amp_s, win, 'omitnan');
med = median(amp_sm(qcFlag), 'omitnan');
mad_ = median(abs(amp_sm(qcFlag) - med), 'omitnan');
thresh = med + 2 * mad_;

[pks, locs] = findpeaks(amp_sm, 'MinPeakHeight', thresh, 'MinPeakDistance', win);

waveEvents = struct('frame', {}, 'time_min', {}, 'phi_shape', {}, ...
                    'phi_vel', {}, 'theta_nuc', {}, ...
                    'dPhi_shape', {}, 'dPhi_vel', {}, 'amp', {}, ...
                    'primary', {});
for ii = 1:numel(locs)
    fr = locs(ii);
    waveEvents(ii).frame      = fr;
    waveEvents(ii).time_min   = time_min(fr);
    waveEvents(ii).phi_shape  = phi_shape(fr);
    waveEvents(ii).phi_vel    = phi_vel(fr);
    waveEvents(ii).theta_nuc  = theta_nuc(fr);
    waveEvents(ii).dPhi_shape = dPhi_shape(fr);
    waveEvents(ii).dPhi_vel   = dPhi_vel(fr);
    waveEvents(ii).amp        = pks(ii);
    waveEvents(ii).primary    = primary;
end
fprintf('Detected %d discrete wave events using %s amplitude (threshold=%.3g).\n', ...
        numel(waveEvents), primary, thresh);

%% -------- CIRCULAR STATS ON dPhi --------
% Prefer shape-axis dPhi; fall back to velocity-axis dPhi when shape is absent.
if strcmp(primary, 'shape')
    dPhi_for_stats = arrayfun(@(e) e.dPhi_shape, waveEvents);
else
    dPhi_for_stats = arrayfun(@(e) e.dPhi_vel,   waveEvents);
end
validE = isfinite(dPhi_for_stats);
if any(validE)
    dPhiE = dPhi_for_stats(validE);
    mean_resultant = mean(exp(1i * dPhiE));
    circ_mean = angle(mean_resultant);
    circ_R    = abs(mean_resultant);
    % Rayleigh test: p ~ exp(-n*R^2) for large n
    n = numel(dPhiE);
    rayleigh_p = exp(-n * circ_R^2) * (1 + (2*n*circ_R^2 - (n*circ_R^2)^2) / (4*n));
    rayleigh_p = min(max(rayleigh_p, 0), 1);
    fprintf('Event-level dPhi_%s: circ mean = %.1f deg, R = %.3f, Rayleigh p ~ %.3g (n=%d)\n', ...
        primary, rad2deg(circ_mean), circ_R, rayleigh_p, n);
else
    circ_mean  = NaN;
    circ_R     = NaN;
    rayleigh_p = NaN;
end

%% -------- SAVE ENRICHED RESULTS --------
fprintf('Saving enriched results to %s ...\n', resultsFile);
save(resultsFile, 'nucleusXY', 'phi_shape', 'amp_shape', 'kStar_shape', ...
     'phi_vel', 'amp_vel', 'theta_nuc', 'dPhi_shape', 'dPhi_vel', ...
     'kymo_slope', 'waveEvents', 'circ_mean', 'circ_R', 'rayleigh_p', ...
     '-append');

%% -------- PLOTS --------
fprintf('Generating plots in %s ...\n', outDir);
plot_wave_axis(resultsFile, outDir);

fprintf('Done. See %s\n', outDir);

%% ========================= HELPERS =======================================
function y = wrap_pi(x)
% Wrap angle to [-pi, pi]
y = mod(x + pi, 2*pi) - pi;
end

function y = wrap_pi_half(x)
% Wrap angle to [-pi/2, pi/2] (for undirected axes)
y = mod(x + pi/2, pi) - pi/2;
end
