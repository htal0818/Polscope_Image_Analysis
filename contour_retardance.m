% contour_retardance.m
% Measure retardance (nm) at the oocyte contour from LC Polscope image stacks.
%
% This script:
%   1. Loads LC Polscope retardance image stacks (multipage TIFF, folder of
%      TIFFs, or 4-state raw Polscope channels)
%   2. Converts raw pixel values to retardance in nm using the Polscope
%      retardance ceiling
%   3. Detects the oocyte outer boundary using heavy Gaussian blur + Otsu
%      thresholding + morphological cleanup + bwboundaries
%   4. Samples retardance (nm) along the detected contour at each time point
%   5. Computes angle-averaged radial retardance profiles from the center
%      outward through the cortex
%   6. Generates:
%        - Retardance kymograph (angle vs time at the cortex)
%        - Radial profiles (center -> cortex -> outside)
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
px_per_um = 6.25;      % pixels per micron (adjust for your objective/camera)

% --- Retardance calibration ---
retardance_ceiling_nm = 50;   % Polscope retardance ceiling (nm)
bit_depth = 16;               % image bit depth (16-bit = 0..65535)

% --- Boundary detection parameters ---
% Uses radial gradient peak method: for each angle from center, find the
% steepest intensity drop (cortex-to-background edge). Robust to internal
% oocyte contrast variations.
sigmaGrad    = 5;       % Gaussian blur sigma (px) before gradient (noise reduction)
nRays        = 360;     % number of radial rays for boundary detection
minEdgeFrac  = 0.3;     % search for edge in outer portion of ray (0.3 = outer 70%)

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
maxPixVal = 2^bit_depth - 1;  % 65535 for 16-bit

% Time axis
time_sec = (0:nFrames-1)' * dt_sec;
time_min = time_sec / 60;

% Angular bins for kymograph
angles_deg = linspace(0, 360, nThetaBins);
thetaBinEdges = linspace(0, 2*pi, nThetaBins+1);

% Radial axis: generous upper bound, trimmed after processing
maxRadius_um = 150;
radialAxis_um = 0 : radialStep_um : maxRadius_um;
nRadial = numel(radialAxis_um);
radialAxis_px = radialAxis_um / um_per_px;

%% ========================== PREALLOCATE ===================================
% Kymograph: retardance at boundary vs angle over time
kymo = nan(nFrames, nThetaBins);

% Contour statistics per frame (in nm)
contourMean = nan(nFrames, 1);
contourStd  = nan(nFrames, 1);
contourMax  = nan(nFrames, 1);
contourMin  = nan(nFrames, 1);

% Radial profiles: [nFrames x nRadial] (in nm)
radialProfiles = nan(nFrames, nRadial);

% Boundary tracking
centroidXY = nan(nFrames, 2);
meanRadius_px = nan(nFrames, 1);


%% ========================== MAIN LOOP =====================================
fprintf('Processing %d frames...\n', nFrames);
tic;

for fr = 1:nFrames

    %% ----- Load image and convert to retardance (nm) -----
    Iraw = readFrame(fr);
    if doCrop
        Iraw = imcrop(Iraw, cropRect);
    end
    [H, W] = size(Iraw);

    % Convert raw pixel values to retardance in nm
    Iret = (Iraw / maxPixVal) * retardance_ceiling_nm;

    %% ----- Boundary detection (radial peak intensity method) -----
    % Strategy: the cortex is the brightest ring in retardance images.
    % For each radial ray from center outward, find the peak intensity
    % in the outer portion — that's directly ON the cortex.

    I_blur = imgaussfilt(Iraw, sigmaGrad);

    % Rough center estimate: centroid of above-background pixels
    I_norm = I_blur / max(I_blur(:));
    roughMask = I_norm > 0.15;  % very permissive, just for centroid
    roughMask = imfill(roughMask, 'holes');
    props = regionprops(roughMask, 'Area', 'Centroid');
    if ~isempty(props)
        [~, iMax] = max([props.Area]);
        roughCenter = props(iMax).Centroid;  % [x, y]
    else
        roughCenter = [W/2, H/2];
    end
    cx0 = roughCenter(1);
    cy0 = roughCenter(2);

    % Maximum radial extent to search (stay within image)
    maxR_px = floor(min([cx0-1, W-cx0, cy0-1, H-cy0])) - 1;
    if maxR_px < 20
        maxR_px = round(min(H, W) / 2) - 1;
    end

    % Cast radial rays and find peak INTENSITY in outer portion of each
    rayAngles = linspace(0, 2*pi, nRays+1);
    rayAngles(end) = [];
    rSamples = 1:maxR_px;
    edgeR = nan(nRays, 1);

    Fraw = griddedInterpolant({1:H, 1:W}, I_blur, 'linear', 'nearest');

    for ri = 1:nRays
        ang = rayAngles(ri);
        xs = cx0 + rSamples * cos(ang);
        ys = cy0 + rSamples * sin(ang);

        % Clip to image bounds
        inBounds = xs >= 1 & xs <= W & ys >= 1 & ys <= H;
        if sum(inBounds) < 10; continue; end

        ivals = Fraw(ys(inBounds), xs(inBounds));
        rvals = rSamples(inBounds);

        % Search for peak intensity in the outer portion of the ray
        % (skip inner region to avoid internal bright spots)
        startIdx = max(1, round(numel(ivals) * minEdgeFrac));
        [~, iPeak] = max(ivals(startIdx:end));
        iPeak = iPeak + startIdx - 1;
        edgeR(ri) = rvals(iPeak);
    end

    % Convert edge points to Cartesian
    validRays = ~isnan(edgeR);
    xb = cx0 + edgeR(validRays) .* cos(rayAngles(validRays)');
    yb = cy0 + edgeR(validRays) .* sin(rayAngles(validRays)');

    if numel(xb) < 20
        fprintf('  Frame %d: too few edge points (%d), skipping.\n', fr, numel(xb));
        continue;
    end

    %% ----- Circle fit for center & radius -----
    [R_fit, xc, yc] = circfit(xb, yb);
    centroidXY(fr,:) = [xc, yc];
    meanRadius_px(fr) = R_fit;

    %% ----- Retardance along the boundary (kymograph row) -----
    % Angle of each boundary point relative to center
    th = atan2(yb - yc, xb - xc);
    th(th < 0) = th(th < 0) + 2*pi;

    % Interpolate retardance (nm) at boundary pixel locations
    F = griddedInterpolant({1:H, 1:W}, Iret, 'linear', 'nearest');
    ib = F(yb, xb);  % retardance in nm at each boundary point

    % Bin by angle
    bin = discretize(th, thetaBinEdges);
    valid = ~isnan(bin);

    row = nan(1, nThetaBins);
    if any(valid)
        sumBins = accumarray(bin(valid), ib(valid), [nThetaBins 1], @nanmean, NaN);
        row = sumBins';
    end

    % Fill missing bins via circular interpolation
    bad = isnan(row);
    if any(bad) && sum(~bad) >= 2
        goodIdx = find(~bad);
        xw = [goodIdx - nThetaBins, goodIdx, goodIdx + nThetaBins];
        vw = [row(goodIdx), row(goodIdx), row(goodIdx)];
        row(bad) = interp1(xw, vw, find(bad), 'linear', 'extrap');
    end

    kymo(fr,:) = row;

    %% ----- Contour retardance statistics (nm) -----
    contourMean(fr) = mean(ib, 'omitnan');
    contourStd(fr)  = std(ib, 'omitnan');
    contourMax(fr)  = max(ib);
    contourMin(fr)  = min(ib);

    %% ----- Angle-averaged radial profile (center outward, in nm) -----
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
        imagesc(Iret); colormap gray; axis image; hold on;

        % Show detected edge points
        plot(xb, yb, 'r.', 'MarkerSize', 4);

        % Show circle fit
        theta_circ = linspace(0, 2*pi, 200);
        plot(xc + R_fit*cos(theta_circ), yc + R_fit*sin(theta_circ), ...
            'c-', 'LineWidth', 1.5);
        plot(xc, yc, 'g+', 'MarkerSize', 12, 'LineWidth', 2);

        title(sprintf('Frame %d / %d  (R=%.0f px = %.0f um)  —  %s', ...
            fr, nFrames, R_fit, R_fit*um_per_px, frameName(fr)), ...
            'Interpreter', 'none');
        cb = colorbar; cb.Label.String = 'Retardance (nm)';
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
cb.Label.String = 'Retardance (nm)';
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
if ~isnan(medianR_um)
    xline(medianR_um, 'r--', 'Cortex', 'LineWidth', 1.5, ...
        'LabelOrientation', 'aligned', 'FontSize', 10);
end
xlabel('Distance from center (\mum)', 'FontSize', 12);
ylabel('Retardance (nm)', 'FontSize', 12);
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
plot(meanRadius_px * um_per_px, time_min, 'r-', 'LineWidth', 1.5);
xlabel('Distance from center (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Retardance: Radial Distance vs Time', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (nm)';
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
ylabel('Retardance (nm)', 'FontSize', 11);
title('Mean Contour Retardance Over Time', 'FontSize', 13);
grid on;
legend('Mean', '\pm1 SD', 'Location', 'best');

subplot(2,1,2);
plot(time_min, contourMax, 'r-', 'LineWidth', 1.2); hold on;
plot(time_min, contourMin, 'b-', 'LineWidth', 1.2);
plot(time_min, contourMean, 'k-', 'LineWidth', 1.5);
xlabel('Time (min)', 'FontSize', 11);
ylabel('Retardance (nm)', 'FontSize', 11);
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
results.kymo                  = kymo;
results.angles_deg            = angles_deg;
results.radialProfiles_nm     = radialProfiles_trim;
results.radialAxis_um         = radialAxis_um_trim;
results.contourMean_nm        = contourMean;
results.contourStd_nm         = contourStd;
results.contourMax_nm         = contourMax;
results.contourMin_nm         = contourMin;
results.centroidXY            = centroidXY;
results.meanRadius_px         = meanRadius_px;
results.meanRadius_um         = meanRadius_px * um_per_px;
results.time_sec              = time_sec;
results.time_min              = time_min;
results.nFrames               = nFrames;
results.dt_sec                = dt_sec;
results.px_per_um             = px_per_um;
results.retardance_ceiling_nm = retardance_ceiling_nm;
results.bit_depth             = bit_depth;
results.nThetaBins            = nThetaBins;
results.nAngleSamples         = nAngleSamples;
results.inputMode             = inputMode;
results.sigmaGrad             = sigmaGrad;
results.nRays                 = nRays;
results.minEdgeFrac           = minEdgeFrac;

save(fullfile(outDir, 'contour_retardance_results.mat'), '-struct', 'results');
fprintf('Saved results to: %s\n', fullfile(outDir, 'contour_retardance_results.mat'));

%% ========================== SUMMARY =======================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('Frames processed: %d\n', sum(~isnan(contourMean)));
fprintf('Duration: %.1f min\n', max(time_min));
fprintf('Median oocyte radius: %.1f um (%.0f px)\n', medianR_um, nanmedian(meanRadius_px));
fprintf('Retardance ceiling: %.0f nm (%d-bit)\n', retardance_ceiling_nm, bit_depth);
fprintf('Mean contour retardance: %.2f +/- %.2f nm\n', ...
    mean(contourMean, 'omitnan'), std(contourMean, 'omitnan'));
fprintf('Outputs saved to: %s\n', outDir);
fprintf('=============================\n');
