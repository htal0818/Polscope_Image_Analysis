% contour_retardance.m
% Measure retardance at the oocyte contour and generate radial profiles over time.
%
% This script:
%   1. Loads LC Polscope retardance image stacks (multipage TIFF, folder of
%      TIFFs, or 4-state raw Polscope channels)
%   2. Detects the oocyte boundary using curvature-based polar polynomial
%      fitting (same methodology as boundary_flows.m)
%   3. Samples retardance along the detected contour at each time point
%   4. Computes angle-averaged radial intensity profiles centered on the
%      boundary (cortex) at each frame
%   5. Generates:
%        - Radial profiles over time (distance from cortex vs retardance)
%        - Mean contour retardance time series
%        - Overlay visualizations of detected boundary
%
% Handles most LC Polscope image stacks:
%   - Single-channel retardance TIFFs (pre-computed by Polscope software)
%   - 4-state raw Polscope data (State1..State4 combined)
%   - Multipage TIFF stacks or folders of individual TIFFs
%
% REQUIREMENTS:
%   - Image Processing Toolbox
%   - circfit.m (included in this repository)

clear all; close all; clc

%% ========================== USER INPUTS ==================================
% --- INPUT MODE ---
% Set inputMode to one of:
%   'retardance'   : pre-computed retardance images (single-channel TIFFs)
%   'four_state'   : raw 4-state Polscope data (State1..State4)
%   'multipage'    : single multipage TIFF file
inputMode = 'retardance';

% --- For 'retardance' mode: folder of single-channel retardance TIFFs ---
base_dir = '/path/to/your/data/Pos0/';
retardance_pattern = '*1_Retardance*';

% --- For 'four_state' mode: folder with State1..State4 images ---
% (uses base_dir above)
% state_patterns = {'*State1*', '*State2*', '*State3*', '*State4*'};

% --- For 'multipage' mode: path to a single multipage TIFF ---
% multipage_path = '/path/to/stack.tif';

% --- Optional crop (set doCrop=false to use full image) ---
doCrop = false;
cropRect = [500 500 2000 2000];  % [x y w h] in pixels

% --- Timing & calibration ---
dt_sec    = 15;        % seconds per frame
px_per_um = 6.25;      % pixels per micron

% --- Boundary detection parameters ---
sigmaBlur  = 1.0;      % Gaussian blur sigma (px) for boundary detection
threshFrac = 0.85;     % mask threshold: I < threshFrac * mean(I)
polyOrder  = 50;        % polynomial order for polar r(theta) fit
nBoundary  = 500;       % number of evenly spaced boundary samples

% --- Quality control ---
minAreaFrac  = 0.05;    % reject if mask area < fraction of image area
maxEccentric = 0.95;    % reject if too eccentric
minSolidity  = 0.85;    % reject if solidity too low

% --- Radial profile parameters ---
radialRange_um = [-5, 10];   % [inside, outside] of boundary in microns
                              % negative = inside cell, positive = outside
radialStep_um  = 0.25;       % radial sampling step (microns)
nAngleSamples  = 360;        % angular resolution for angle-averaged profile

% --- Output ---
outDir = fullfile(base_dir, 'contour_retardance_out');

% --- Visualization ---
saveOverlays     = true;       % save boundary overlay images
overlayEveryN    = 10;         % save overlay every N frames
profileEveryN    = 20;         % save individual radial profile every N frames

% ============================================================================
%% ========================== LOAD IMAGE LIST ================================

switch inputMode
    case 'retardance'
        d = dir(fullfile(base_dir, retardance_pattern));
        if isempty(d)
            error('No retardance images found matching "%s" in %s', ...
                retardance_pattern, base_dir);
        end
        [~, sortIdx] = sort({d.name});
        d = d(sortIdx);
        nFrames = numel(d);
        readFrame = @(t) double(imread(fullfile(d(t).folder, d(t).name)));
        frameName = @(t) d(t).name;

    case 'four_state'
        state_patterns = {'*State1*', '*State2*', '*State3*', '*State4*'};
        ds = cell(1,4);
        for si = 1:4
            ds{si} = dir(fullfile(base_dir, state_patterns{si}));
            [~, idx] = sort({ds{si}.name});
            ds{si} = ds{si}(idx);
        end
        nFrames = min(cellfun(@numel, ds));
        readFrame = @(t) double(imread(fullfile(ds{1}(t).folder, ds{1}(t).name))) ...
                       + double(imread(fullfile(ds{2}(t).folder, ds{2}(t).name))) ...
                       + double(imread(fullfile(ds{3}(t).folder, ds{3}(t).name))) ...
                       + double(imread(fullfile(ds{4}(t).folder, ds{4}(t).name)));
        frameName = @(t) ds{1}(t).name;

    case 'multipage'
        info = imfinfo(multipage_path);
        nFrames = numel(info);
        readFrame = @(t) double(imread(multipage_path, t));
        frameName = @(t) sprintf('page %d', t);

    otherwise
        error('Unknown inputMode: %s. Use ''retardance'', ''four_state'', or ''multipage''.', inputMode);
end

fprintf('Found %d frames (mode: %s)\n', nFrames, inputMode);

%% ========================== SETUP =========================================
if ~exist(outDir, 'dir'); mkdir(outDir); end
if saveOverlays
    overlayDir = fullfile(outDir, 'overlays');
    if ~exist(overlayDir, 'dir'); mkdir(overlayDir); end
end

um_per_px = 1 / px_per_um;

% Radial axis (in microns, centered on boundary)
radialAxis_um = radialRange_um(1) : radialStep_um : radialRange_um(2);
nRadial = numel(radialAxis_um);
radialAxis_px = radialAxis_um / um_per_px;

% Angles for radial sampling
sampleAngles = linspace(0, 2*pi, nAngleSamples+1);
sampleAngles(end) = [];

% Time axis
time_sec = (0:nFrames-1)' * dt_sec;
time_min = time_sec / 60;

%% ========================== PREALLOCATE ===================================
% Radial profiles: [nFrames x nRadial]
radialProfiles = nan(nFrames, nRadial);

% Contour statistics per frame
contourMean = nan(nFrames, 1);
contourStd  = nan(nFrames, 1);
contourMax  = nan(nFrames, 1);
contourMin  = nan(nFrames, 1);

% Boundary tracking
centroidXY = nan(nFrames, 2);
meanRadius = nan(nFrames, 1);

%% ========================== MAIN LOOP =====================================
fprintf('Processing %d frames...\n', nFrames);
tic;

for fr = 1:nFrames

    %% ----- Load image -----
    Iraw = readFrame(fr);
    if doCrop
        Iraw = imcrop(Iraw, cropRect);
    end
    [H, W] = size(Iraw);
    I = mat2gray(Iraw);  % normalized for boundary detection

    %% ----- Boundary detection (curvature-based polar polynomial) -----

    % Step 1: Coarse threshold mask
    Ib = imgaussfilt(I, sigmaBlur);
    BW = Ib < threshFrac * mean2(Ib);
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, round(minAreaFrac * H * W));

    % Morphological cleanup
    se = strel('diamond', 5);
    BW = imclose(BW, se);
    BW = imfill(BW, 'holes');

    % Keep largest connected component
    L = bwlabel(BW, 8);
    if max(L(:)) >= 1
        S = regionprops(L, 'Area', 'Centroid', 'Eccentricity', 'Solidity');
        [~, iMax] = max([S.Area]);
        BW = (L == iMax);

        % Quality control
        if S(iMax).Eccentricity > maxEccentric || S(iMax).Solidity < minSolidity
            fprintf('  Frame %d: QC warning (ecc=%.2f, sol=%.2f)\n', ...
                fr, S(iMax).Eccentricity, S(iMax).Solidity);
        end
    else
        % Fallback: gradient-based detection
        [Gmag, ~] = imgradient(Ib);
        thrG = max(2*mean(Gmag(:)), prctile(Gmag(:), 85));
        BW = imfill(bwareaopen(Gmag >= thrG, round(minAreaFrac*H*W)), 'holes');
        L = bwlabel(BW, 8);
        if max(L(:)) >= 1
            S = regionprops(L, 'Area', 'Centroid');
            [~, iMax] = max([S.Area]);
            BW = (L == iMax);
        else
            fprintf('  Frame %d: no boundary found, skipping.\n', fr);
            continue;
        end
    end

    % Step 2: Edge detection & circle fit for center
    edgeBW = bwperim(BW);
    [ey, ex] = find(edgeBW);

    if numel(ex) < 10
        fprintf('  Frame %d: too few edge points (%d), skipping.\n', fr, numel(ex));
        continue;
    end

    [R_fit, xc, yc] = circfit(ex, ey);
    centroidXY(fr,:) = [xc, yc];
    meanRadius(fr) = R_fit;

    % Step 3: Polar coordinate transformation of edge points
    theta_edge = atan2(ey - yc, ex - xc);
    theta_edge(theta_edge < 0) = theta_edge(theta_edge < 0) + 2*pi;
    r_edge = sqrt((ex - xc).^2 + (ey - yc).^2);

    % Step 4: Two-pass polynomial fit for r(theta) with wrap-around
    % First pass: standard polynomial fit
    [theta_sort, si] = sort(theta_edge);
    r_sort = r_edge(si);

    % Remove outliers (points too far from circle fit)
    rResid = abs(r_sort - R_fit);
    rThresh = 2 * median(rResid);
    keep = rResid < rThresh;
    theta_clean = theta_sort(keep);
    r_clean = r_sort(keep);

    if numel(r_clean) < polyOrder + 1
        polyOrd = max(5, numel(r_clean) - 1);
    else
        polyOrd = polyOrder;
    end

    % Wrap-around: extend data cyclically for smooth fit at 0/2pi
    theta_ext = [theta_clean - 2*pi; theta_clean; theta_clean + 2*pi];
    r_ext = [r_clean; r_clean; r_clean];
    p = polyfit(theta_ext, r_ext, polyOrd);

    % Evaluate smooth boundary
    theta_smooth = linspace(0, 2*pi, nBoundary+1);
    theta_smooth(end) = [];
    r_smooth = polyval(p, theta_smooth);

    % Enforce minimum radius (prevent self-intersection)
    r_smooth = max(r_smooth, R_fit * 0.3);

    % Convert back to Cartesian
    xBnd = xc + r_smooth .* cos(theta_smooth);
    yBnd = yc + r_smooth .* sin(theta_smooth);

    %% ----- Sample retardance along radial lines from boundary -----
    % For each angle, cast a radial line from center through boundary
    % and sample retardance at distances relative to the boundary crossing

    % Interpolant for retardance values
    F = griddedInterpolant({1:H, 1:W}, Iraw, 'linear', 'nearest');

    % Angle-averaged radial profile
    profile_sum   = zeros(1, nRadial);
    profile_count = zeros(1, nRadial);

    for ai = 1:nAngleSamples
        ang = sampleAngles(ai);

        % Boundary radius at this angle
        r_bnd = polyval(p, ang);
        r_bnd = max(r_bnd, R_fit * 0.3);

        % Sample points along radial direction, centered on boundary
        % positive radialAxis_px = outward (away from center)
        % negative radialAxis_px = inward (toward center)
        r_samples = r_bnd + radialAxis_px;

        % Convert to image coordinates
        xs = xc + r_samples * cos(ang);
        ys = yc + r_samples * sin(ang);

        % Clip to image bounds
        valid = xs >= 1 & xs <= W & ys >= 1 & ys <= H;

        if any(valid)
            vals = F(ys(valid), xs(valid));
            profile_sum(valid) = profile_sum(valid) + vals(:)';
            profile_count(valid) = profile_count(valid) + 1;
        end
    end

    % Compute angle-averaged profile
    validBins = profile_count > 0;
    radialProfiles(fr, validBins) = profile_sum(validBins) ./ profile_count(validBins);

    % Contour statistics (at radial distance = 0, i.e., the boundary itself)
    [~, zeroBin] = min(abs(radialAxis_um));
    contourMean(fr) = radialProfiles(fr, zeroBin);

    % Also compute stats from all boundary sample points
    xBndSample = xc + r_smooth .* cos(theta_smooth);
    yBndSample = yc + r_smooth .* sin(theta_smooth);
    xBndSample = max(1, min(W, xBndSample));
    yBndSample = max(1, min(H, yBndSample));
    bndVals = F(yBndSample(:), xBndSample(:));
    contourStd(fr) = std(bndVals, 'omitnan');
    contourMax(fr) = max(bndVals);
    contourMin(fr) = min(bndVals);

    %% ----- Save overlay -----
    if saveOverlays && (mod(fr, overlayEveryN) == 1 || fr == 1)
        fig = figure('Visible', 'off', 'Position', [100 100 800 600]);
        imagesc(Iraw); colormap gray; axis image; hold on;
        plot([xBnd, xBnd(1)], [yBnd, yBnd(1)], 'r-', 'LineWidth', 1.5);
        plot(xc, yc, 'g+', 'MarkerSize', 12, 'LineWidth', 2);
        title(sprintf('Frame %d / %d  —  %s', fr, nFrames, frameName(fr)), ...
            'Interpreter', 'none');
        colorbar;
        exportgraphics(gca, fullfile(overlayDir, sprintf('overlay_%04d.png', fr)));
        close(fig);
    end

    %% ----- Save individual radial profile -----
    if mod(fr, profileEveryN) == 1 || fr == 1
        fig = figure('Visible', 'off', 'Position', [100 100 700 400]);
        plot(radialAxis_um, radialProfiles(fr,:), 'b-', 'LineWidth', 1.5);
        hold on;
        xline(0, 'r--', 'Cortex', 'LineWidth', 1.2, 'LabelOrientation', 'aligned');
        xlabel('Distance from cortex (\mum)');
        ylabel('Retardance (a.u.)');
        title(sprintf('Radial profile — Frame %d (t = %.0f s)', fr, time_sec(fr)));
        grid on;
        exportgraphics(gca, fullfile(outDir, sprintf('radial_profile_%04d.png', fr)));
        close(fig);
    end

    % Progress
    if mod(fr, 25) == 0 || fr == nFrames
        fprintf('  Processed %d / %d frames (%.0f%%)\n', fr, nFrames, 100*fr/nFrames);
    end
end

elapsed = toc;
fprintf('Done! %.1f sec total (%.2f sec/frame)\n\n', elapsed, elapsed/nFrames);

%% ========================== GENERATE PLOTS ================================
fprintf('Generating summary plots...\n');

% --- Plot 1: Radial profiles over time (waterfall / stacked) ---
% Select evenly spaced frames to display
nDisplay = min(20, nFrames);
displayFrames = unique(round(linspace(1, nFrames, nDisplay)));

fig1 = figure('Position', [100 100 900 600]);
cmap = parula(numel(displayFrames));
hold on;
for i = 1:numel(displayFrames)
    fr = displayFrames(i);
    if all(isnan(radialProfiles(fr,:))); continue; end
    plot(radialAxis_um, radialProfiles(fr,:), '-', 'Color', cmap(i,:), ...
        'LineWidth', 1.2, 'DisplayName', sprintf('t=%.0fs', time_sec(fr)));
end
xline(0, 'k--', 'LineWidth', 1.5);
xlabel('Distance from cortex (\mum)', 'FontSize', 12);
ylabel('Retardance (a.u.)', 'FontSize', 12);
title('Radial Retardance Profiles Over Time', 'FontSize', 14);
cb = colorbar;
colormap(gca, parula(numel(displayFrames)));
clim([time_sec(displayFrames(1)), time_sec(displayFrames(end))]);
cb.Label.String = 'Time (s)';
cb.Label.FontSize = 11;
grid on;
legend('show', 'Location', 'eastoutside', 'FontSize', 7);
exportgraphics(fig1, fullfile(outDir, 'radial_profiles_over_time.png'), 'Resolution', 200);
savefig(fig1, fullfile(outDir, 'radial_profiles_over_time.fig'));
close(fig1);

% --- Plot 2: Radial profile heatmap (distance vs time) ---
fig2 = figure('Position', [100 100 900 500]);
imagesc(radialAxis_um, time_min, radialProfiles);
set(gca, 'YDir', 'normal');
hold on;
xline(0, 'r-', 'LineWidth', 1.5);
xlabel('Distance from cortex (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Retardance: Radial Distance vs Time', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (a.u.)';
cb.Label.FontSize = 11;
exportgraphics(fig2, fullfile(outDir, 'radial_heatmap_distance_vs_time.png'), 'Resolution', 200);
savefig(fig2, fullfile(outDir, 'radial_heatmap_distance_vs_time.fig'));
close(fig2);

% --- Plot 3: Mean contour retardance over time ---
fig3 = figure('Position', [100 100 800 450]);
subplot(2,1,1);
plot(time_min, contourMean, 'b-', 'LineWidth', 1.5);
hold on;
fill([time_min; flipud(time_min)], ...
     [contourMean - contourStd; flipud(contourMean + contourStd)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (min)', 'FontSize', 11);
ylabel('Retardance (a.u.)', 'FontSize', 11);
title('Mean Contour Retardance Over Time', 'FontSize', 13);
grid on;
legend('Mean', '\pm1 SD', 'Location', 'best');

subplot(2,1,2);
plot(time_min, contourStd, 'r-', 'LineWidth', 1.5);
xlabel('Time (min)', 'FontSize', 11);
ylabel('Std Dev (a.u.)', 'FontSize', 11);
title('Contour Retardance Variability', 'FontSize', 13);
grid on;

exportgraphics(fig3, fullfile(outDir, 'contour_retardance_timeseries.png'), 'Resolution', 200);
savefig(fig3, fullfile(outDir, 'contour_retardance_timeseries.fig'));
close(fig3);

% --- Plot 4: Boundary radius over time ---
fig4 = figure('Position', [100 100 700 350]);
plot(time_min, meanRadius * um_per_px, 'k-', 'LineWidth', 1.5);
xlabel('Time (min)', 'FontSize', 11);
ylabel('Mean radius (\mum)', 'FontSize', 11);
title('Oocyte Radius Over Time', 'FontSize', 13);
grid on;
exportgraphics(fig4, fullfile(outDir, 'oocyte_radius_over_time.png'), 'Resolution', 200);
close(fig4);

%% ========================== SAVE DATA =====================================
results = struct();
results.radialProfiles  = radialProfiles;
results.radialAxis_um   = radialAxis_um;
results.contourMean     = contourMean;
results.contourStd      = contourStd;
results.contourMax      = contourMax;
results.contourMin      = contourMin;
results.centroidXY      = centroidXY;
results.meanRadius_px   = meanRadius;
results.time_sec        = time_sec;
results.time_min        = time_min;
results.nFrames         = nFrames;
results.dt_sec          = dt_sec;
results.px_per_um       = px_per_um;
results.radialRange_um  = radialRange_um;
results.radialStep_um   = radialStep_um;
results.nAngleSamples   = nAngleSamples;
results.inputMode       = inputMode;
results.sigmaBlur       = sigmaBlur;
results.threshFrac      = threshFrac;
results.polyOrder       = polyOrder;
results.nBoundary       = nBoundary;

save(fullfile(outDir, 'contour_retardance_results.mat'), '-struct', 'results');
fprintf('Saved results to: %s\n', fullfile(outDir, 'contour_retardance_results.mat'));

%% ========================== SUMMARY =======================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('Frames processed: %d\n', sum(~isnan(contourMean)));
fprintf('Duration: %.1f min\n', max(time_min));
fprintf('Radial range: [%.1f, %.1f] um from cortex\n', radialRange_um(1), radialRange_um(2));
fprintf('Mean contour retardance: %.2f +/- %.2f\n', ...
    mean(contourMean, 'omitnan'), std(contourMean, 'omitnan'));
fprintf('Outputs saved to: %s\n', outDir);
fprintf('=============================\n');
