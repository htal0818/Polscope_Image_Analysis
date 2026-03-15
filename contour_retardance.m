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
% Heavy blur washes out internal oocyte structure so Otsu finds the gross
% egg shape. Only used for mask creation — all measurements use raw data.
sigmaBlur   = 20;      % Gaussian blur sigma (px) for segmentation mask
closeRadius = 25;      % morphological close disk radius (px)
minArea     = 5000;    % minimum object area (px^2) to reject debris
boundaryInset_px = 10; % shift boundary inward (px) onto cortical ring center

% --- Threshold mode ---
% 'otsu'       : automatic Otsu threshold on blurred image (default)
% 'fixed'      : fixed intensity threshold on raw image
% 'percentile' : threshold at a percentile of raw image intensity
thresholdMode       = 'otsu';
fixedThreshold      = 500;     % raw pixel value for 'fixed' mode
percentileThreshold = 30;      % percentile for 'percentile' mode (pixels ABOVE this)

% --- Radial profile parameters ---
radialStep_um  = 0.5;  % radial sampling step (microns)
nAngleSamples  = 360;  % angular resolution for radial profiles

% --- Outside-in radial profile parameters ---
nFourier       = 25;   % number of Fourier harmonics for boundary fit
nBoundaryPts   = 500;  % number of uniformly spaced boundary points
maxDepth_um    = 50;   % how far inward from cortex to sample (microns)
depthStep_um   = 0.5;  % step size along inward normals (microns)

% --- Mask caching (reuse mask between frames for speed) ---
useMaskCaching          = true;   % master switch: false = recalculate every frame
cacheIntensityThreshold = 0.02;   % reuse mask if mean intensity change < 2%
cacheForceRecalcEveryN  = 25;     % force full recalc every N frames (drift correction)

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

% Radial axis (center-out, kept for legacy compatibility)
maxRadius_um = 150;
radialAxis_um = 0 : radialStep_um : maxRadius_um;
nRadial = numel(radialAxis_um);
radialAxis_px = radialAxis_um / um_per_px;

% Depth axis (outside-in: 0 = cortex, increasing = deeper into oocyte)
depthAxis_um = 0 : depthStep_um : maxDepth_um;
nDepth = numel(depthAxis_um);
depthAxis_px = depthAxis_um / um_per_px;

%% ========================== PREALLOCATE ===================================
% Kymograph: retardance at boundary vs angle over time
kymo = nan(nFrames, nThetaBins);

% Contour statistics per frame (in nm)
contourMean = nan(nFrames, 1);
contourStd  = nan(nFrames, 1);
contourMax  = nan(nFrames, 1);
contourMin  = nan(nFrames, 1);

% Radial profiles: [nFrames x nRadial] (in nm) — center-out (legacy)
radialProfiles = nan(nFrames, nRadial);

% Outside-in profiles: normal-based [nFrames x nDepth] (angle-averaged)
normalProfiles = nan(nFrames, nDepth);

% Outside-in profiles: distance transform [nFrames x nDepth] (angle-averaged)
distProfiles = nan(nFrames, nDepth);

% Outside-in 2D map: normal-based [nBoundaryPts x nDepth] per frame (last frame stored)
% Full kymograph-style: [nFrames x nDepth] for distance transform
distKymo = nan(nFrames, nDepth);

% Boundary tracking
centroidXY = nan(nFrames, 2);
meanRadius_px = nan(nFrames, 1);

% --- Cache state for smart mask reuse ---
cache_prevBW      = [];    % previous binary mask
cache_prevMeanInt = [];    % previous mean intensity inside mask
cacheHitCount     = 0;     % diagnostic counter
cacheRecalcCount  = 0;     % diagnostic counter

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

    %% ----- Boundary detection (with smart mask caching) -----

    % Decide whether to recalculate or reuse previous mask
    needsRecalc = true;

    if useMaskCaching && ~isempty(cache_prevBW)
        if mod(fr, cacheForceRecalcEveryN) == 1
            needsRecalc = true;   % forced drift correction
        else
            % Check intensity change inside previous mask
            I_norm_check = Iraw / max(Iraw(:));
            meanIntCurrent = mean(I_norm_check(cache_prevBW), 'omitnan');
            intensityChange = abs(meanIntCurrent - cache_prevMeanInt) / (cache_prevMeanInt + eps);

            if intensityChange < cacheIntensityThreshold
                needsRecalc = false;  % cache hit
            end
        end
    end

    if needsRecalc || ~useMaskCaching
        % --- FULL MASK RECALCULATION ---
        switch thresholdMode
            case 'otsu'
                I_blur = imgaussfilt(Iraw, sigmaBlur);
                I_norm = I_blur / max(I_blur(:));
                Totsu = graythresh(I_norm);
                BW = I_norm > Totsu;

            case 'fixed'
                I_blur = imgaussfilt(Iraw, sigmaBlur);
                BW = I_blur > fixedThreshold;

            case 'percentile'
                I_blur = imgaussfilt(Iraw, sigmaBlur);
                pVal = prctile(I_blur(:), percentileThreshold);
                BW = I_blur > pVal;

            otherwise
                error('Unknown thresholdMode: %s', thresholdMode);
        end

        % Aggressive morphological cleanup
        se = strel('disk', closeRadius);
        BW = imclose(BW, se);
        BW = imfill(BW, 'holes');
        BW = bwareaopen(BW, minArea);

        % Fallback: gradient-based if threshold fails
        if ~any(BW(:))
            [Gmag, ~] = imgradient(I_blur);
            thrG = max(2*mean(Gmag(:)), prctile(Gmag(:), 80));
            BW = Gmag >= thrG;
            BW = imclose(BW, se);
            BW = imfill(BW, 'holes');
            BW = bwareaopen(BW, minArea);
        end

        % Keep largest connected component (fall back to cached mask on failure)
        L = bwlabel(BW, 8);
        if max(L(:)) >= 1
            S = regionprops(L, 'Area', 'Centroid');
            [~, iMax] = max([S.Area]);
            BW = (L == iMax);
        elseif ~isempty(cache_prevBW)
            BW = cache_prevBW;
        else
            fprintf('  Frame %d: no boundary found, skipping.\n', fr);
            continue;
        end

        % Update cache
        if useMaskCaching
            cache_prevBW = BW;
            I_norm_cache = Iraw / max(Iraw(:));
            cache_prevMeanInt = mean(I_norm_cache(BW), 'omitnan');
        end
        cacheRecalcCount = cacheRecalcCount + 1;

    else
        % --- CACHE REUSE (skip segmentation) ---
        BW = cache_prevBW;
        cacheHitCount = cacheHitCount + 1;
    end

    %% ----- Extract boundary contour -----
    B = bwboundaries(BW);
    if isempty(B)
        fprintf('  Frame %d: bwboundaries returned empty, skipping.\n', fr);
        continue;
    end
    [~, iLongest] = max(cellfun(@(p) size(p,1), B));
    bnd = B{iLongest};
    yb = bnd(:,1);
    xb = bnd(:,2);

    %% ----- Circle fit for center & radius -----
    [R_fit, xc, yc] = circfit(xb, yb);
    centroidXY(fr,:) = [xc, yc];
    meanRadius_px(fr) = R_fit;

    % Shrink boundary inward onto cortical ring center
    dx = xb - xc;  dy = yb - yc;
    dist = sqrt(dx.^2 + dy.^2);
    shrink = max(dist - boundaryInset_px, 1) ./ dist;
    xb = xc + dx .* shrink;
    yb = yc + dy .* shrink;

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

    %% ----- Fourier boundary fit (outside-in profiling) -----
    % Fourier series is naturally periodic — no wrap-around artifacts
    % and no Runge's phenomenon (unlike high-order polynomial fit).
    bnd_dx = xb - xc;  bnd_dy = yb - yc;
    bnd_r = sqrt(bnd_dx.^2 + bnd_dy.^2);
    bnd_theta = atan2(bnd_dy, bnd_dx);  % range [-pi, pi]

    % Sort by angle for fitting
    [bnd_theta_sort, sIdx] = sort(bnd_theta);
    bnd_r_sort = bnd_r(sIdx);

    % Build design matrix: r(theta) = a0 + sum_n [an*cos(n*theta) + bn*sin(n*theta)]
    nPts_bnd = numel(bnd_theta_sort);
    A_fourier = ones(nPts_bnd, 2*nFourier + 1);
    for ni = 1:nFourier
        A_fourier(:, 2*ni)   = cos(ni * bnd_theta_sort);
        A_fourier(:, 2*ni+1) = sin(ni * bnd_theta_sort);
    end
    fourier_coeffs = A_fourier \ bnd_r_sort;

    % Evaluate on uniform grid
    polyTheta = linspace(-pi, pi, nBoundaryPts+1);
    polyTheta(end) = [];
    A_eval = ones(1, 2*nFourier + 1);
    A_eval = repmat(A_eval, nBoundaryPts, 1);
    A_eval(:,1) = 1;
    for ni = 1:nFourier
        A_eval(:, 2*ni)   = cos(ni * polyTheta(:));
        A_eval(:, 2*ni+1) = sin(ni * polyTheta(:));
    end
    polyR = (A_eval * fourier_coeffs)';

    % Smooth boundary points in Cartesian
    polyX = xc + polyR .* cos(polyTheta);
    polyY = yc + polyR .* sin(polyTheta);

    % ---- Analytic derivative dr/dtheta from Fourier coefficients ----
    drdtheta = zeros(1, nBoundaryPts);
    for ni = 1:nFourier
        an = fourier_coeffs(2*ni);
        bn = fourier_coeffs(2*ni+1);
        drdtheta = drdtheta - ni * an * sin(ni * polyTheta) ...
                             + ni * bn * cos(ni * polyTheta);
    end

    % Tangent vector in Cartesian: d/dtheta [r*cos(theta), r*sin(theta)]
    tx = drdtheta .* cos(polyTheta) - polyR .* sin(polyTheta);
    ty = drdtheta .* sin(polyTheta) + polyR .* cos(polyTheta);
    tn = sqrt(tx.^2 + ty.^2);
    tx = tx ./ tn;  ty = ty ./ tn;

    % Inward normal = rotate tangent 90 deg clockwise (points toward center)
    % Check: normal should point roughly toward center
    nx = ty;   ny = -tx;
    % Flip normals that point outward (dot with center direction)
    toCenter_x = xc - polyX;  toCenter_y = yc - polyY;
    dot_check = nx .* toCenter_x + ny .* toCenter_y;
    nx(dot_check < 0) = -nx(dot_check < 0);
    ny(dot_check < 0) = -ny(dot_check < 0);

    %% ----- Normal-based outside-in radial profiles -----
    normal_sum   = zeros(1, nDepth);
    normal_count = zeros(1, nDepth);

    for bi = 1:nBoundaryPts
        % Sample along inward normal from this boundary point
        xs_n = polyX(bi) + depthAxis_px * nx(bi);
        ys_n = polyY(bi) + depthAxis_px * ny(bi);

        inBounds = xs_n >= 1 & xs_n <= W & ys_n >= 1 & ys_n <= H;
        if any(inBounds)
            vals = F(ys_n(inBounds), xs_n(inBounds));
            normal_sum(inBounds)   = normal_sum(inBounds)   + vals(:)';
            normal_count(inBounds) = normal_count(inBounds) + 1;
        end
    end
    validN = normal_count > 0;
    normalProfiles(fr, validN) = normal_sum(validN) ./ normal_count(validN);

    %% ----- Distance transform outside-in radial profiles -----
    per = bwperim(BW);
    D = bwdist(per) * um_per_px;  % distance from cortex in microns

    % Mask interior only
    D_interior = D;
    D_interior(~BW) = NaN;

    % Bin all interior pixels by their distance from cortex
    depthBinEdges = [depthAxis_um - depthStep_um/2, depthAxis_um(end) + depthStep_um/2];
    depthBins = discretize(D_interior(:), depthBinEdges);
    retVals = Iret(:);
    validD = ~isnan(depthBins);
    if any(validD)
        distProfiles(fr, :) = accumarray(depthBins(validD), retVals(validD), ...
            [nDepth 1], @nanmean, NaN)';
    end

    % Also store for the depth-vs-time kymograph
    distKymo(fr, :) = distProfiles(fr, :);

    %% ----- Save overlay -----
    if saveOverlays && (fr == 1 || mod(fr, overlayEveryN) == 0)
        fig = figure('Visible', 'off', 'Position', [100 100 800 600]);
        imagesc(Iret); colormap gray; axis image; hold on;

        % Show detected boundary (shrunk inward)
        plot(xb, yb, 'r-', 'LineWidth', 1.2);

        % Show polynomial-fit boundary (smooth)
        plot([polyX polyX(1)], [polyY polyY(1)], 'g-', 'LineWidth', 1.5);

        % Show a few inward normals (every 25th point)
        normalVis_px = 15;  % length of normal arrows in pixels
        for vi = 1:25:nBoundaryPts
            plot([polyX(vi), polyX(vi) + normalVis_px*nx(vi)], ...
                 [polyY(vi), polyY(vi) + normalVis_px*ny(vi)], ...
                 'y-', 'LineWidth', 0.8);
        end

        % Show circle fit
        theta_circ = linspace(0, 2*pi, 200);
        plot(xc + R_fit*cos(theta_circ), yc + R_fit*sin(theta_circ), ...
            'c--', 'LineWidth', 0.8);
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

% --- Plot 6: Outside-in radial profiles (normal-based, selected frames) ---
fig6 = figure('Position', [100 100 900 550]);
cmap6 = parula(numel(displayFrames));
hold on;
for i = 1:numel(displayFrames)
    fr = displayFrames(i);
    if all(isnan(normalProfiles(fr,:))); continue; end
    plot(depthAxis_um, normalProfiles(fr,:), '-', ...
        'Color', cmap6(i,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('t=%.0fs', time_sec(fr)));
end
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Retardance (nm)', 'FontSize', 12);
title('Outside-In Radial Profiles (normal-based)', 'FontSize', 14);
grid on;
legend('show', 'Location', 'eastoutside', 'FontSize', 7);
exportgraphics(fig6, fullfile(outDir, 'radial_profiles_normal_outside_in.png'), 'Resolution', 200);
savefig(fig6, fullfile(outDir, 'radial_profiles_normal_outside_in.fig'));
close(fig6);

% --- Plot 7: Outside-in radial profiles (distance transform, selected frames) ---
fig7 = figure('Position', [100 100 900 550]);
cmap7 = parula(numel(displayFrames));
hold on;
for i = 1:numel(displayFrames)
    fr = displayFrames(i);
    if all(isnan(distProfiles(fr,:))); continue; end
    plot(depthAxis_um, distProfiles(fr,:), '-', ...
        'Color', cmap7(i,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('t=%.0fs', time_sec(fr)));
end
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Retardance (nm)', 'FontSize', 12);
title('Outside-In Radial Profiles (distance transform)', 'FontSize', 14);
grid on;
legend('show', 'Location', 'eastoutside', 'FontSize', 7);
exportgraphics(fig7, fullfile(outDir, 'radial_profiles_dist_outside_in.png'), 'Resolution', 200);
savefig(fig7, fullfile(outDir, 'radial_profiles_dist_outside_in.fig'));
close(fig7);

% --- Plot 8: Depth-vs-time heatmap (distance transform) ---
fig8 = figure('Position', [100 100 900 500]);
imagesc(depthAxis_um, time_min, distKymo);
set(gca, 'YDir', 'normal');
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Retardance: Depth from Cortex vs Time (distance transform)', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (nm)';
cb.Label.FontSize = 11;
exportgraphics(fig8, fullfile(outDir, 'depth_heatmap_dist_vs_time.png'), 'Resolution', 200);
savefig(fig8, fullfile(outDir, 'depth_heatmap_dist_vs_time.fig'));
close(fig8);

% --- Plot 9: Comparison — normal-based vs distance transform (time-averaged) ---
fig9 = figure('Position', [100 100 800 450]);
meanNormal = nanmean(normalProfiles, 1);
meanDist   = nanmean(distProfiles, 1);
plot(depthAxis_um, meanNormal, 'b-', 'LineWidth', 2, 'DisplayName', 'Normal-based');
hold on;
plot(depthAxis_um, meanDist, 'r-', 'LineWidth', 2, 'DisplayName', 'Distance transform');
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Retardance (nm)', 'FontSize', 12);
title('Time-Averaged Outside-In Profiles: Normal vs Distance Transform', 'FontSize', 13);
grid on;
legend('show', 'Location', 'best', 'FontSize', 11);
exportgraphics(fig9, fullfile(outDir, 'comparison_normal_vs_dist.png'), 'Resolution', 200);
savefig(fig9, fullfile(outDir, 'comparison_normal_vs_dist.fig'));
close(fig9);

%% ========================== SAVE DATA =====================================
results = struct();
results.kymo                  = kymo;
results.angles_deg            = angles_deg;
results.radialProfiles_nm     = radialProfiles_trim;
results.radialAxis_um         = radialAxis_um_trim;
results.normalProfiles_nm     = normalProfiles;
results.distProfiles_nm       = distProfiles;
results.distKymo_nm           = distKymo;
results.depthAxis_um          = depthAxis_um;
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
results.thresholdMode         = thresholdMode;
results.nFourier              = nFourier;
results.nBoundaryPts          = nBoundaryPts;
results.maxDepth_um           = maxDepth_um;
results.depthStep_um          = depthStep_um;
results.sigmaBlur             = sigmaBlur;
results.closeRadius           = closeRadius;
results.minArea               = minArea;
results.boundaryInset_px      = boundaryInset_px;

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
if useMaskCaching
    totalCached = cacheHitCount + cacheRecalcCount;
    if totalCached > 0
        fprintf('Mask cache hits: %d / %d (%.1f%%)\n', ...
            cacheHitCount, totalCached, 100*cacheHitCount/totalCached);
    end
end
fprintf('Outputs saved to: %s\n', outDir);
fprintf('=============================\n');
