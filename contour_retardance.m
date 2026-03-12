% contour_retardance.m
% Measure retardance at the oocyte contour from LC Polscope image stacks.
%
% This script:
%   1. Loads LC Polscope retardance image stacks (multipage TIFF, folder of
%      TIFFs, or 4-state raw Polscope channels)
%   2. Detects the oocyte boundary using Otsu thresholding + bwboundaries
%      (same proven approach as image_segment_test.m)
%   3. Samples retardance along the detected contour at each time point
%   4. Computes angle-averaged radial intensity profiles from the center
%      outward through the cortex
%   5. Generates:
%        - Retardance kymograph (angle vs time at the cortex)
%        - Radial profiles (center → cortex → outside)
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
sigmaBlur = 2;         % Gaussian blur sigma (px) for segmentation
minArea   = 200;       % minimum object area (px^2) to keep

% --- Cortical band width (for radial sampling near boundary) ---
bandWidth_um = 3;      % sample retardance within this distance of boundary (um)

% --- Radial profile parameters ---
radialStep_um  = 0.5;  % radial sampling step (microns)
nAngleSamples  = 360;  % angular resolution for radial profiles

% --- Angular binning for kymograph ---
nThetaBins = 100;      % number of angular bins around contour

% --- Output ---
outDir = fullfile(base_dir, 'contour_retardance_out');

% --- Visualization ---
saveOverlays     = true;    % save boundary overlay images
overlayEveryN    = 10;      % save overlay every N frames

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
bandWidth_px = bandWidth_um * px_per_um;

% Time axis
time_sec = (0:nFrames-1)' * dt_sec;
time_min = time_sec / 60;

% Angular bins for kymograph
angles_deg = linspace(0, 360, nThetaBins);
thetaBinEdges = linspace(0, 2*pi, nThetaBins+1);

% Radial axis: from center out to ~1.5x expected oocyte radius
% (will be trimmed after first frame establishes actual radius)
maxRadius_um = 120;  % generous upper bound; will auto-trim
radialAxis_um = 0 : radialStep_um : maxRadius_um;
nRadial = numel(radialAxis_um);
radialAxis_px = radialAxis_um / um_per_px;

%% ========================== PREALLOCATE ===================================
% Kymograph: retardance at boundary vs angle over time
kymo = nan(nFrames, nThetaBins);

% Contour statistics per frame
contourMean = nan(nFrames, 1);
contourStd  = nan(nFrames, 1);
contourMax  = nan(nFrames, 1);
contourMin  = nan(nFrames, 1);

% Radial profiles: [nFrames x nRadial]
radialProfiles = nan(nFrames, nRadial);

% Boundary tracking
centroidXY = nan(nFrames, 2);
meanRadius_px = nan(nFrames, 1);

% For fallback: previous frame's mask
prevBW = [];
prevC  = [];

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
    I = im2double(Iraw / max(Iraw(:)));  % normalize to [0,1]

    %% ----- Boundary detection (Otsu + bwboundaries) -----
    % Follows the proven approach from image_segment_test.m:
    % oocyte is BRIGHT in retardance images, so threshold picks bright regions

    I_blur = imgaussfilt(I, sigmaBlur);

    % Otsu threshold
    Totsu = graythresh(I_blur);
    BW = I_blur > Totsu;
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, minArea);

    % Fallback: gradient-based if Otsu fails
    if ~any(BW(:))
        [Gmag, ~] = imgradient(I_blur);
        thrG = max(2*mean(Gmag(:)), prctile(Gmag(:), 80));
        BW = Gmag >= thrG;
        BW = imfill(BW, 'holes');
        BW = bwareaopen(BW, minArea);
    end

    % Keep largest connected component
    L = bwlabel(BW, 8);
    if max(L(:)) >= 1
        S = regionprops(L, 'Area', 'Centroid');
        [~, iMax] = max([S.Area]);
        BW = (L == iMax);
        C = S(iMax).Centroid;  % [x, y]
    elseif ~isempty(prevBW)
        % Reuse previous frame's mask
        BW = prevBW;
        C = prevC;
    else
        fprintf('  Frame %d: no boundary found, skipping.\n', fr);
        continue;
    end

    prevBW = BW;
    prevC = C;
    centroidXY(fr,:) = C;

    %% ----- Extract boundary contour -----
    B = bwboundaries(BW);
    if isempty(B)
        fprintf('  Frame %d: bwboundaries returned empty, skipping.\n', fr);
        continue;
    end

    % Pick the longest boundary (largest object)
    [~, iLongest] = max(cellfun(@(p) size(p,1), B));
    bnd = B{iLongest};  % [row, col] = [y, x]
    yb = bnd(:,1);
    xb = bnd(:,2);

    %% ----- Circle fit for center & radius -----
    [R_fit, xc, yc] = circfit(xb, yb);
    centroidXY(fr,:) = [xc, yc];
    meanRadius_px(fr) = R_fit;

    %% ----- Retardance along the boundary (kymograph row) -----
    % Compute angle of each boundary point relative to center
    th = atan2(yb - yc, xb - xc);
    th(th < 0) = th(th < 0) + 2*pi;

    % Interpolate retardance at boundary pixel locations
    F = griddedInterpolant({1:H, 1:W}, Iraw, 'linear', 'nearest');
    ib = F(yb, xb);  % retardance intensity at each boundary point

    % Bin by angle
    bin = discretize(th, thetaBinEdges);
    valid = ~isnan(bin);

    row = nan(1, nThetaBins);
    if any(valid)
        sumBins   = accumarray(bin(valid), ib(valid), [nThetaBins 1], @nanmean, NaN);
        row = sumBins';
    end

    % Fill missing bins via circular interpolation
    bad = isnan(row);
    if any(bad) && sum(~bad) >= 2
        goodIdx = find(~bad);
        % Wrap for circular interpolation
        xw = [goodIdx - nThetaBins, goodIdx, goodIdx + nThetaBins];
        vw = [row(goodIdx), row(goodIdx), row(goodIdx)];
        row(bad) = interp1(xw, vw, find(bad), 'linear', 'extrap');
    end

    kymo(fr,:) = row;

    %% ----- Contour retardance statistics -----
    contourMean(fr) = mean(ib, 'omitnan');
    contourStd(fr)  = std(ib, 'omitnan');
    contourMax(fr)  = max(ib);
    contourMin(fr)  = min(ib);

    %% ----- Angle-averaged radial profile (center outward) -----
    sampleAngles = linspace(0, 2*pi, nAngleSamples+1);
    sampleAngles(end) = [];

    profile_sum   = zeros(1, nRadial);
    profile_count = zeros(1, nRadial);

    for ai = 1:nAngleSamples
        ang = sampleAngles(ai);

        % Sample along ray from center outward
        xs = xc + radialAxis_px * cos(ang);
        ys = yc + radialAxis_px * sin(ang);

        % Keep points inside image
        inBounds = xs >= 1 & xs <= W & ys >= 1 & ys <= H;

        if any(inBounds)
            vals = F(ys(inBounds), xs(inBounds));
            profile_sum(inBounds)   = profile_sum(inBounds)   + vals(:)';
            profile_count(inBounds) = profile_count(inBounds) + 1;
        end
    end

    validR = profile_count > 0;
    radialProfiles(fr, validR) = profile_sum(validR) ./ profile_count(validR);

    %% ----- Save overlay -----
    if saveOverlays && (fr == 1 || mod(fr, overlayEveryN) == 0)
        fig = figure('Visible', 'off', 'Position', [100 100 800 600]);
        imagesc(Iraw); colormap gray; axis image; hold on;
        plot(xb, yb, 'r-', 'LineWidth', 1.2);
        plot(xc, yc, 'g+', 'MarkerSize', 12, 'LineWidth', 2);

        % Show cortical band
        theta_circ = linspace(0, 2*pi, 200);
        plot(xc + R_fit*cos(theta_circ), yc + R_fit*sin(theta_circ), ...
            'c--', 'LineWidth', 0.8);
        title(sprintf('Frame %d / %d  —  %s', fr, nFrames, frameName(fr)), ...
            'Interpreter', 'none');
        colorbar;
        exportgraphics(gca, fullfile(overlayDir, sprintf('overlay_%04d.png', fr)));
        close(fig);
    end

    % Progress
    if mod(fr, 25) == 0 || fr == nFrames
        fprintf('  Processed %d / %d frames (%.0f%%)\n', fr, nFrames, 100*fr/nFrames);
    end
end

elapsed = toc;
fprintf('Done! %.1f sec total (%.2f sec/frame)\n\n', elapsed, elapsed/nFrames);

%% ========================== TRIM RADIAL AXIS ==============================
% Trim radial profiles to 1.5x the median oocyte radius
medianR_um = nanmedian(meanRadius_px) * um_per_px;
trimIdx = find(radialAxis_um <= medianR_um * 1.5, 1, 'last');
if isempty(trimIdx); trimIdx = nRadial; end
radialAxis_um_trim = radialAxis_um(1:trimIdx);
radialProfiles_trim = radialProfiles(:, 1:trimIdx);

%% ========================== GENERATE PLOTS ================================
fprintf('Generating summary plots...\n');

% --- Plot 1: Retardance kymograph (angle vs time at cortex) ---
fig1 = figure('Position', [100 100 900 500]);
imagesc(angles_deg, time_min, kymo);
set(gca, 'YDir', 'normal');
xlabel('Angle around cortex (deg)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Retardance at Cortex (angle vs time)', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (a.u.)';
cb.Label.FontSize = 11;
exportgraphics(fig1, fullfile(outDir, 'kymograph_retardance_vs_angle.png'), 'Resolution', 200);
savefig(fig1, fullfile(outDir, 'kymograph_retardance_vs_angle.fig'));
close(fig1);

% --- Plot 2: Radial profiles over time (selected frames) ---
nDisplay = min(15, nFrames);
displayFrames = unique(round(linspace(1, nFrames, nDisplay)));

fig2 = figure('Position', [100 100 900 550]);
cmap = parula(numel(displayFrames));
hold on;
for i = 1:numel(displayFrames)
    fr = displayFrames(i);
    if all(isnan(radialProfiles_trim(fr,:))); continue; end
    plot(radialAxis_um_trim, radialProfiles_trim(fr,:), '-', ...
        'Color', cmap(i,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('t=%.0fs', time_sec(fr)));
end
% Mark median cortex position
if ~isnan(medianR_um)
    xline(medianR_um, 'r--', 'Cortex', 'LineWidth', 1.5, ...
        'LabelOrientation', 'aligned', 'FontSize', 10);
end
xlabel('Distance from center (\mum)', 'FontSize', 12);
ylabel('Retardance (a.u.)', 'FontSize', 12);
title('Radial Retardance Profiles (center \rightarrow cortex)', 'FontSize', 14);
grid on;
legend('show', 'Location', 'eastoutside', 'FontSize', 7);
exportgraphics(fig2, fullfile(outDir, 'radial_profiles_over_time.png'), 'Resolution', 200);
savefig(fig2, fullfile(outDir, 'radial_profiles_over_time.fig'));
close(fig2);

% --- Plot 3: Radial profile heatmap (distance vs time) ---
fig3 = figure('Position', [100 100 900 500]);
imagesc(radialAxis_um_trim, time_min, radialProfiles_trim);
set(gca, 'YDir', 'normal');
hold on;
% Mark cortex position per frame
plot(meanRadius_px * um_per_px, time_min, 'r-', 'LineWidth', 1.5);
xlabel('Distance from center (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Retardance: Radial Distance vs Time', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (a.u.)';
cb.Label.FontSize = 11;
legend('Cortex boundary', 'Location', 'northeast');
exportgraphics(fig3, fullfile(outDir, 'radial_heatmap_distance_vs_time.png'), 'Resolution', 200);
savefig(fig3, fullfile(outDir, 'radial_heatmap_distance_vs_time.fig'));
close(fig3);

% --- Plot 4: Mean contour retardance over time ---
fig4 = figure('Position', [100 100 800 500]);

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
plot(time_min, contourMax, 'r-', 'LineWidth', 1.2); hold on;
plot(time_min, contourMin, 'b-', 'LineWidth', 1.2);
plot(time_min, contourMean, 'k-', 'LineWidth', 1.5);
xlabel('Time (min)', 'FontSize', 11);
ylabel('Retardance (a.u.)', 'FontSize', 11);
title('Contour Retardance Range', 'FontSize', 13);
legend('Max', 'Min', 'Mean', 'Location', 'best');
grid on;

exportgraphics(fig4, fullfile(outDir, 'contour_retardance_timeseries.png'), 'Resolution', 200);
savefig(fig4, fullfile(outDir, 'contour_retardance_timeseries.fig'));
close(fig4);

% --- Plot 5: Oocyte radius over time ---
fig5 = figure('Position', [100 100 700 350]);
plot(time_min, meanRadius_px * um_per_px, 'k-', 'LineWidth', 1.5);
xlabel('Time (min)', 'FontSize', 11);
ylabel('Mean radius (\mum)', 'FontSize', 11);
title('Oocyte Radius Over Time', 'FontSize', 13);
grid on;
exportgraphics(fig5, fullfile(outDir, 'oocyte_radius_over_time.png'), 'Resolution', 200);
close(fig5);

%% ========================== SAVE DATA =====================================
results = struct();
results.kymo            = kymo;
results.angles_deg      = angles_deg;
results.radialProfiles  = radialProfiles_trim;
results.radialAxis_um   = radialAxis_um_trim;
results.contourMean     = contourMean;
results.contourStd      = contourStd;
results.contourMax      = contourMax;
results.contourMin      = contourMin;
results.centroidXY      = centroidXY;
results.meanRadius_px   = meanRadius_px;
results.meanRadius_um   = meanRadius_px * um_per_px;
results.time_sec        = time_sec;
results.time_min        = time_min;
results.nFrames         = nFrames;
results.dt_sec          = dt_sec;
results.px_per_um       = px_per_um;
results.nThetaBins      = nThetaBins;
results.nAngleSamples   = nAngleSamples;
results.inputMode       = inputMode;
results.sigmaBlur       = sigmaBlur;
results.minArea         = minArea;

save(fullfile(outDir, 'contour_retardance_results.mat'), '-struct', 'results');
fprintf('Saved results to: %s\n', fullfile(outDir, 'contour_retardance_results.mat'));

%% ========================== SUMMARY =======================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('Frames processed: %d\n', sum(~isnan(contourMean)));
fprintf('Duration: %.1f min\n', max(time_min));
fprintf('Median oocyte radius: %.1f um\n', medianR_um);
fprintf('Mean contour retardance: %.2f +/- %.2f\n', ...
    mean(contourMean, 'omitnan'), std(contourMean, 'omitnan'));
fprintf('Outputs saved to: %s\n', outDir);
fprintf('=============================\n');
