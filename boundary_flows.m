% Mesure tangential surface flows from PIV data 


% building off the main idea from SCW_flows, this script will calculate the
% tangential vector field, or all vectors perpendicular to the oocyte
% cortex. the magniutude and directionality (CW or CCW) detail the strength
% of cortical contractions. 




%% SCW_tangential_kymograph_curvature.m
% Tangential cortical flow kymograph from PIVlab + strict oocyte boundary via curvature.
% Uses curvature-based boundary reconstruction methodology adapted from kymograph.m
% as implemented in SCW_flows_curvature.m for PIV data.
%
% REQUIREMENTS:
% - Image Processing Toolbox (imfill, bwdist, etc.)
% - circfit.m for circle fitting to boundary points

clear all; close all; clc

%% ========================== USER INPUTS ==================================
% --- PolScope image folders / patterns (like your existing code) ---
base_dir = '/Users/hridaytalreja/Desktop/Jan_data_2026/jan_20_2026_FSW_and_eggs_50msexp_15sint_20x_50nmceiling/eggs/SMS_2026_0120_1518_1/Pos0/';

% If you have 4 states per frame (best for boundary):
s1 = fullfile(base_dir,'*State1*');
s2 = fullfile(base_dir,'*State2*');
s3 = fullfile(base_dir,'*State3*');
s4 = fullfile(base_dir,'*State4*');

% If you only have State1, set useFourStates=false and fill s1 accordingly.
useFourStates = true;

d1 = dir(s1);
if useFourStates
    d2 = dir(s2); d3 = dir(s3); d4 = dir(s4);
end

% --- PIVlab output .mat file ---
pivMatFile = '/Users/hridaytalreja/Desktop/Jan_data_2026/jan_20_2026_FSW_and_eggs_50msexp_15sint_20x_50nmceiling/eggs/SMS_2026_0120_1518_1/Pos0/Jan21_2026_PIV/PIVlab_output.mat';

% --- Optional crop (match your workflows) ---
doCrop = true;
cropRect = [500 500 2000 2000];  % [x y w h] in pixels (like your script)

% --- Units / timing ---
px_per_um = 6.25 /2;     % px/um
dt_sec    = 15;       % sec/frame

% --- PIV velocity units ---
pivVelUnit = 'px_per_frame';  % 'px_per_frame' or 'px_per_sec'

% --- Cortical band definition (STRICT): pixels INSIDE boundary ---
% Example: bandOuter=1, bandInner=20  => sample points 1–20 px inside cortex
% Wider band captures more cortical flow signal
bandOuterPx = 1;
bandInnerPx = 20;

% --- Curvature-based boundary parameters (from kymograph.m methodology) ---
sigmaBlur       = 1.0;           % pre-blur (pixels); keep small to preserve edge
threshFrac      = 0.85;          % mask = I < threshFrac*mean2(I)
se              = strel('diamond',5);
polyOrder       = 50;            % polyfit order for r(theta) (matches kymograph)
nBoundary       = 500;           % samples along boundary for smooth polygon

% --- Theta binning (like your kymograph) ---
nThetaBins = 101;                 % like your theta=linspace(0,2*pi,101)
thetaEdges = linspace(0,2*pi,nThetaBins);  % bin centers style
thetaCenters = thetaEdges;        % treat as centers for output
thetaBinEdges = linspace(0,2*pi,nThetaBins+1);

% --- Boundary smoothing / spline sampling ---
nDenseSpline = 5000;   % dense samples along boundary for smooth nearest-point lookup (increased density)

% --- Quality control thresholds ---
minAreaFrac  = 0.05;   % reject if mask area < frac of image area
maxEccentric = 0.95;   % reject if region is too eccentric (bad segmentation)
minSolidity  = 0.85;   % reject if region solidity too low (noisy boundary)

% ============================================================================
% --- SMART CACHING FOR MASK REUSABILITY ---
% Set useMaskCaching = false to disable and recalculate every frame
% ============================================================================
useMaskCaching = true;   % MASTER SWITCH: set to false to disable all caching

% Caching strategy (only active if useMaskCaching = true):
cacheIntensityThreshold = 0.02;  % Reuse mask if mean intensity change < 2%
cacheForceRecalcEveryN  = 25;    % Force full recalc every N frames (drift correction)
cacheUseHintCenter      = true;  % Use previous center as circle fit initialization

% --- Outputs ---
outDir = fullfile(base_dir, 'tangential_kymo_out');
if ~exist(outDir,'dir'); mkdir(outDir); end

% Create subdirectories for organized outputs
outDir_kymographs = fullfile(outDir, 'kymographs');
outDir_quiver = fullfile(outDir, 'quiver_overlays');
outDir_snapshots = fullfile(outDir, 'snapshots');
outDir_qc = fullfile(outDir, 'qc');

if ~exist(outDir_kymographs,'dir'); mkdir(outDir_kymographs); end
if ~exist(outDir_quiver,'dir'); mkdir(outDir_quiver); end
if ~exist(outDir_snapshots,'dir'); mkdir(outDir_snapshots); end
if ~exist(outDir_qc,'dir'); mkdir(outDir_qc); end

saveEveryN_QC = 20;

% ============================================================================
% --- VISUALIZATION SETTINGS ---
% ============================================================================
makeQuiverOverlays = true;   % Create quiver plot overlays on oocyte boundary
quiverEveryNFrames = 10;     % Save quiver overlay every N frames
quiverSubsample = 3;         % Subsample theta bins for quiver (every Nth bin)
quiverScale = 1.5;           % Arrow length scaling factor

makeEnhancedKymographs = true;  % Create magnitude + signed directional kymographs
makeSnapshotPlots = true;       % Create detailed snapshot visualizations
snapshotFrames = [];            % Specific frames to visualize (empty = auto-select)
snapshotEveryNFrames = 2;       % Save tangential flow analysis every N frames

%% ========================== LOAD PIV DATA =================================
S = load(pivMatFile);
[Xc, Yc, Uc, Vc] = pickPIVFields(S);

nFramesPIV = numel(Uc);
nFramesImg = numel(d1);
nFrames = min(nFramesPIV, nFramesImg);
fprintf('Frames: images=%d, piv=%d, using=%d\n', nFramesImg, nFramesPIV, nFrames);

%% ========================== PREALLOCATE ===================================
Vtheta_kymo = nan(nFrames, nThetaBins);   % mean tangential velocity per theta bin (um/s)
Npts_kymo   = zeros(nFrames, nThetaBins); % counts per bin
centroidXY  = nan(nFrames,2);
areaMask    = nan(nFrames,1);
qcFlag      = false(nFrames,1);           % true if frame passes QC
RADIUS_OF_CURVATURE = nan(nFrames, numel(thetaCenters)-1);  % radius of curvature per theta bin

% --- Cache variables (for smart mask reusability) ---
cache_prevBW = [];           % Previous mask
cache_prevCenter = [];       % Previous [xc, yc]
cache_prevRadius = [];       % Previous mean radius
cache_prevMeanInt = [];      % Previous mean intensity in mask region
cache_prevPoly = [];         % Previous boundary polygon
cache_prevRADIUS = [];       % Previous curvature array
cacheHitCount = 0;           % Diagnostic: number of cache reuses
cacheRecalcCount = 0;        % Diagnostic: number of full recalculations

% --- Visualization data storage ---
visPolySeq = cell(nFrames,1);     % Store polygon for each frame (for quiver plots)
visImageSeq = cell(nFrames,1);    % Store normalized images (for overlays)

%% ========================== MAIN LOOP =====================================
for fr = 1:nFrames

    %% ----- Load PolScope intensity image used for boundary -----
    a1 = double(imread(fullfile(d1(fr).folder, d1(fr).name)));
    if useFourStates
        a2 = double(imread(fullfile(d2(fr).folder, d2(fr).name)));
        a3 = double(imread(fullfile(d3(fr).folder, d3(fr).name)));
        a4 = double(imread(fullfile(d4(fr).folder, d4(fr).name)));
        Iraw = a1 + a2 + a3 + a4;
    else
        Iraw = a1;
    end

    if doCrop
        Iraw = imcrop(Iraw, cropRect);
    end

    I = mat2gray(Iraw);

    %% =========================================================================
    %% SMART CACHING BLOCK - Set useMaskCaching = false at top to disable
    %% =========================================================================

    % --- Decide if we need full recalculation or can reuse previous mask ---
    needsRecalc = true;  % default: recalculate
    cacheReasonStr = 'first frame';

    if useMaskCaching && ~isempty(cache_prevBW)
        % Check cache validity criteria

        % Criterion 1: Force recalculation every N frames (drift correction)
        if mod(fr, cacheForceRecalcEveryN) == 1
            needsRecalc = true;
            cacheReasonStr = sprintf('forced recalc (every %d frames)', cacheForceRecalcEveryN);
        else
            % Criterion 2: Check intensity change in previous mask region
            meanIntCurrent = mean(I(cache_prevBW), 'omitnan');
            intensityChange = abs(meanIntCurrent - cache_prevMeanInt) / (cache_prevMeanInt + eps);

            if intensityChange < cacheIntensityThreshold
                needsRecalc = false;  % CACHE HIT!
                cacheReasonStr = sprintf('cache hit (Δintensity=%.3f%%)', intensityChange*100);
            else
                needsRecalc = true;  % CACHE MISS
                cacheReasonStr = sprintf('intensity changed %.3f%%', intensityChange*100);
            end
        end
    end

    % --- Execute mask calculation (with or without caching) ---
    if needsRecalc || ~useMaskCaching
        % FULL RECALCULATION
        if ~useMaskCaching
            cacheReasonStr = 'caching disabled';
        end

        % Optionally use previous center as initialization hint
        centerHint = [];
        if useMaskCaching && cacheUseHintCenter && ~isempty(cache_prevCenter)
            centerHint = cache_prevCenter;
        end

        [BW, stats, poly, RADIUS, xc, yc] = make_oocyte_mask_curvature(I, ...
            sigmaBlur, threshFrac, se, polyOrder, nBoundary, thetaCenters, ...
            minAreaFrac, maxEccentric, minSolidity, centerHint);

        if isempty(stats)
            fprintf('Frame %d: mask failed QC (%s).\n', fr, cacheReasonStr);
            continue;
        end

        % Update cache
        if useMaskCaching
            cache_prevBW = BW;
            cache_prevCenter = [xc, yc];
            cache_prevRadius = mean(RADIUS, 'omitnan');
            cache_prevMeanInt = mean(I(BW), 'omitnan');
            cache_prevPoly = poly;
            cache_prevRADIUS = RADIUS;
        end

        cacheRecalcCount = cacheRecalcCount + 1;

    else
        % CACHE REUSE
        BW = cache_prevBW;
        poly = cache_prevPoly;
        RADIUS = cache_prevRADIUS;
        xc = cache_prevCenter(1);
        yc = cache_prevCenter(2);

        % Quick QC check on cached mask
        stats = regionprops(BW, 'Area', 'Eccentricity', 'Centroid', 'Solidity');
        if isempty(stats)
            fprintf('Frame %d: cached mask invalid, forcing recalc.\n', fr);
            % Force recalculation next iteration
            cache_prevBW = [];
            continue;
        end
        [~, ii] = max([stats.Area]);
        stats = stats(ii);

        cacheHitCount = cacheHitCount + 1;
    end

    %% =========================================================================
    %% END SMART CACHING BLOCK
    %% =========================================================================

    % Store results (common path for both cached and recalculated)
    if isempty(stats)
        fprintf('Frame %d: mask failed QC.\n', fr);
        continue;
    end

    qcFlag(fr)    = true;
    centroidXY(fr,:) = [xc, yc];
    areaMask(fr)  = stats.Area;
    RADIUS_OF_CURVATURE(fr,:) = RADIUS;

    %% ----- Get PIV frame -----
    X = Xc{fr}; Y = Yc{fr}; U = Uc{fr}; V = Vc{fr};
    if doCrop
        % If you cropped images, you MUST also shift PIV coordinates accordingly.
        % PIVlab X,Y are in image coordinates. Cropping moves origin by cropRect(1:2).
        X = X - cropRect(1);
        Y = Y - cropRect(2);
    end

    %% ----- Compute tangential v_theta(theta) using curvature-based boundary -----
    [vBins, nBinsCount, dbgFlow] = tangential_from_boundary_curvature( ...
        X, Y, U, V, BW, poly, [xc, yc], ...
        bandOuterPx, bandInnerPx, nThetaBins, thetaBinEdges, ...
        px_per_um, dt_sec, pivVelUnit, nDenseSpline);

    Vtheta_kymo(fr,:) = vBins;
    Npts_kymo(fr,:)   = nBinsCount;

    %% ----- Store visualization data -----
    visPolySeq{fr} = poly;
    visImageSeq{fr} = I;

    %% ----- QC overlays (saved periodically) -----
    if mod(fr, saveEveryN_QC) == 1
        qcFig = figure('Visible','off'); imshow(I,[]); hold on;
        plot(dbgFlow.xDense, dbgFlow.yDense, 'LineWidth', 2);
        scatter(dbgFlow.sampleX, dbgFlow.sampleY, 8, 'filled');
        title(sprintf('Frame %d boundary + sampled PIV band points', fr));
        exportgraphics(qcFig, fullfile(outDir_qc, sprintf('QC_boundary_band_fr%04d.png', fr)), 'Resolution', 250);
        close(qcFig);
    end

    if mod(fr,50)==0
        fprintf('Processed frame %d/%d (%s)\n', fr, nFrames, cacheReasonStr);
    end
end

%% ========================== CACHE DIAGNOSTICS ================================
if useMaskCaching
    totalProcessed = cacheHitCount + cacheRecalcCount;
    cacheHitRate = 100 * cacheHitCount / max(1, totalProcessed);
    fprintf('\n=== SMART CACHING PERFORMANCE ===\n');
    fprintf('Cache hits:        %d/%d (%.1f%%)\n', cacheHitCount, totalProcessed, cacheHitRate);
    fprintf('Recalculations:    %d/%d (%.1f%%)\n', cacheRecalcCount, totalProcessed, 100-cacheHitRate);
    fprintf('Estimated speedup: %.1fx\n', 1 + (cacheHitCount * 9) / max(1, totalProcessed));
    fprintf('=================================\n\n');
else
    fprintf('\nSmart caching was DISABLED (useMaskCaching=false)\n\n');
end

%% ========================== VISUALIZATIONS ===================================
time_min = (0:nFrames-1) * (dt_sec/60);

fprintf('Generating visualizations...\n');

% --- 1) BASIC SIGNED KYMOGRAPH (original) ---
fig1 = figure('Position', [100 100 1000 600]);
imagesc(1:nThetaBins, time_min, Vtheta_kymo);
axis tight;
xlabel('Angular Bin Number', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Tangential Cortical Flow v_\theta(\theta,t) - Signed (μm/s)', 'FontSize', 14);
colormap(gca, 'parula');
cb = colorbar;
ylabel(cb, 'v_\theta (μm/s)', 'FontSize', 11);
set(gca, 'FontSize', 11);
exportgraphics(fig1, fullfile(outDir_kymographs,'kymograph_vtheta_signed.png'), 'Resolution', 300);

if makeEnhancedKymographs
    % --- 2) MAGNITUDE KYMOGRAPH (absolute values) ---
    fig2 = figure('Position', [100 100 1000 600]);
    imagesc(1:nThetaBins, time_min, abs(Vtheta_kymo));
    axis tight;
    xlabel('Angular Bin Number', 'FontSize', 12);
    ylabel('Time (min)', 'FontSize', 12);
    title('Tangential Flow Magnitude |v_\theta(\theta,t)| (μm/s)', 'FontSize', 14);
    colormap(gca, 'hot');
    cb = colorbar;
    ylabel(cb, '|v_\theta| (μm/s)', 'FontSize', 11);
    set(gca, 'FontSize', 11);
    exportgraphics(fig2, fullfile(outDir_kymographs,'kymograph_vtheta_magnitude.png'), 'Resolution', 300);

    % --- 3) DIRECTIONAL KYMOGRAPH (diverging colormap) ---
    fig3 = figure('Position', [100 100 1000 600]);
    imagesc(1:nThetaBins, time_min, Vtheta_kymo);
    axis tight;
    xlabel('Angular Bin Number', 'FontSize', 12);
    ylabel('Time (min)', 'FontSize', 12);
    title('Tangential Flow Directionality (μm/s)', 'FontSize', 14);

    % Diverging colormap: blue (negative/clockwise) to red (positive/counterclockwise)
    colormap(gca, redblue(256));

    % Symmetric color limits around zero
    vmax = max(abs(Vtheta_kymo(:)), [], 'omitnan');
    if ~isnan(vmax) && vmax > 0
        clim([-vmax, vmax]);
    end

    cb = colorbar;
    ylabel(cb, 'v_\theta (μm/s): red=CCW, blue=CW', 'FontSize', 10);
    set(gca, 'FontSize', 11);
    exportgraphics(fig3, fullfile(outDir_kymographs,'kymograph_vtheta_directional.png'), 'Resolution', 300);

    fprintf('  - Saved 3 kymograph variants\n');
end

% --- 4) QUIVER OVERLAYS (tangential flow vectors on boundary) ---
if makeQuiverOverlays
    fprintf('  - Generating quiver overlays...\n');

    % Compute global velocity range for consistent colorbar across all frames
    allVtheta = Vtheta_kymo(qcFlag, :);
    vmin_global = min(allVtheta(:), [], 'omitnan');
    vmax_global = max(allVtheta(:), [], 'omitnan');
    vabs_max = max(abs([vmin_global, vmax_global]));

    % Use symmetric color limits for diverging colormap
    clim_quiver = [-vabs_max, vabs_max];

    % Create diverging colormap (blue=CW/negative, white=zero, red=CCW/positive)
    cmap_diverge = redblue(256);

    for fr = 1:quiverEveryNFrames:nFrames
        if ~qcFlag(fr), continue; end

        I = visImageSeq{fr};
        poly = visPolySeq{fr};
        if isempty(I) || isempty(poly), continue; end

        % Get center and tangential velocities
        xc = centroidXY(fr, 1);
        yc = centroidXY(fr, 2);
        vtheta = Vtheta_kymo(fr, :);

        % Create figure
        figQ = figure('Visible', 'off', 'Position', [100 100 900 800]);
        imshow(I, []); hold on; axis on;

        % Plot boundary (thinner, cyan for contrast)
        plot(poly(:,1), poly(:,2), 'c-', 'LineWidth', 1);

        % Compute quiver positions and vectors
        thetaSubsample = 1:quiverSubsample:nThetaBins;
        nQuiver = numel(thetaSubsample);

        xq = zeros(nQuiver, 1);
        yq = zeros(nQuiver, 1);
        uq = zeros(nQuiver, 1);
        vq = zeros(nQuiver, 1);
        vq_color = zeros(nQuiver, 1);  % velocity values for coloring

        for k = 1:nQuiver
            idx = thetaSubsample(k);
            theta_k = thetaCenters(idx);
            vtheta_k = vtheta(idx);

            if isnan(vtheta_k)
                vq_color(k) = NaN;
                continue;
            end

            % Position on boundary (approximate using polar coordinates from center)
            R_mean = mean(RADIUS_OF_CURVATURE(fr, :), 'omitnan');
            if isnan(R_mean), R_mean = 100; end

            xq(k) = xc + R_mean * cos(theta_k);
            yq(k) = yc + R_mean * sin(theta_k);

            % Tangent vector (perpendicular to radial direction)
            tx = -sin(theta_k);
            ty = cos(theta_k);

            % Scale by velocity magnitude (always positive for arrow length)
            % Direction is encoded by the tangent vector sign
            arrow_scale = quiverScale * abs(vtheta_k) / px_per_um;  % pixels

            % Flip direction for negative (CW) velocities
            if vtheta_k < 0
                tx = -tx;
                ty = -ty;
            end

            uq(k) = arrow_scale * tx;
            vq(k) = arrow_scale * ty;
            vq_color(k) = vtheta_k;  % store actual velocity for color
        end

        % Remove NaN entries
        valid = ~isnan(vq_color) & (xq ~= 0 | yq ~= 0);
        xq = xq(valid); yq = yq(valid);
        uq = uq(valid); vq = vq(valid);
        vq_color = vq_color(valid);

        % Draw color-coded arrows individually
        % Map velocity to colormap index
        nColors = size(cmap_diverge, 1);
        colorIdx = round((vq_color - clim_quiver(1)) / (clim_quiver(2) - clim_quiver(1)) * (nColors - 1)) + 1;
        colorIdx = max(1, min(nColors, colorIdx));  % clamp to valid range

        for k = 1:numel(xq)
            arrowColor = cmap_diverge(colorIdx(k), :);

            % Draw arrow as line with arrowhead
            quiver(xq(k), yq(k), uq(k), vq(k), 0, ...
                'Color', arrowColor, 'LineWidth', 1.5, 'MaxHeadSize', 2);
        end

        % Add scatter plot to drive the colorbar (invisible dots)
        hScatter = scatter(xq, yq, 1, vq_color, 'filled');
        set(hScatter, 'MarkerFaceAlpha', 0);  % make invisible

        % Set colormap and color limits
        colormap(gca, cmap_diverge);
        caxis(clim_quiver);

        % Add colorbar
        cb = colorbar;
        ylabel(cb, 'v_\theta (\mum/s): red=CCW, blue=CW', 'FontSize', 10, 'Color', 'w');
        set(cb, 'Color', 'w');  % white tick labels

        title(sprintf('Frame %d: Tangential Flow (t=%.2f min)', fr, time_min(fr)), ...
            'FontSize', 12, 'Color', 'w');

        exportgraphics(figQ, fullfile(outDir_quiver, sprintf('quiver_tangential_fr%04d.png', fr)), ...
            'Resolution', 200);
        close(figQ);
    end

    fprintf('  - Saved %d quiver overlays to %s\n', numel(1:quiverEveryNFrames:nFrames), outDir_quiver);
end

% --- 5) SNAPSHOT DETAILED VISUALIZATIONS (Tangential Flow Analysis) ---
if makeSnapshotPlots
    fprintf('  - Generating tangential flow analysis plots (every %d frames)...\n', snapshotEveryNFrames);

    for fr = 1:snapshotEveryNFrames:nFrames
        if ~qcFlag(fr), continue; end

        I = visImageSeq{fr};
        poly = visPolySeq{fr};
        if isempty(I) || isempty(poly), continue; end

        xc = centroidXY(fr, 1);
        yc = centroidXY(fr, 2);
        vtheta = Vtheta_kymo(fr, :);
        radius_curv = RADIUS_OF_CURVATURE(fr, :);

        % Create multi-panel figure
        figSnap = figure('Position', [50 50 1400 800]);

        % Panel 1: Image with boundary and quiver
        subplot(2,3,1);
        imshow(I, []); hold on;
        plot(poly(:,1), poly(:,2), 'y-', 'LineWidth', 2);
        plot(xc, yc, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
        title(sprintf('Frame %d (t=%.2f min)', fr, time_min(fr)));

        % Panel 2: Tangential velocity profile
        subplot(2,3,2);
        plot(1:nThetaBins, vtheta, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
        xlabel('Angular Bin Number'); ylabel('v_\theta (μm/s)');
        title('Tangential Velocity Profile');
        grid on;
        xlim([1 nThetaBins]);

        % Panel 3: Curvature profile
        subplot(2,3,3);
        plot(1:nThetaBins-1, radius_curv, 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8);
        xlabel('Angular Bin Number'); ylabel('Radius of Curvature (px)');
        title('Boundary Curvature');
        grid on;
        xlim([1 nThetaBins-1]);

        % Panel 4: Polar plot of tangential velocity
        subplot(2,3,4);
        polarplot(thetaCenters, abs(vtheta), 'b-', 'LineWidth', 2);
        title('|v_\theta| Magnitude (Polar)');

        % Panel 5: Signed polar plot
        subplot(2,3,5);
        % For signed polar plot, show as two colors
        pos_idx = vtheta >= 0;
        neg_idx = vtheta < 0;
        polarplot(thetaCenters(pos_idx), vtheta(pos_idx), 'r.', 'MarkerSize', 8); hold on;
        polarplot(thetaCenters(neg_idx), abs(vtheta(neg_idx)), 'b.', 'MarkerSize', 8);
        title('v_\theta Directionality');
        legend({'CCW (+)', 'CW (-)'}, 'Location', 'northeast');  % Top right corner

        % Panel 6: Statistics
        subplot(2,3,6); axis off;
        vtheta_valid = vtheta(~isnan(vtheta));
        if ~isempty(vtheta_valid)
            stats_text = {
                sprintf('Frame: %d', fr)
                sprintf('Time: %.2f min', time_min(fr))
                ''
                'Tangential Flow Statistics:'
                sprintf('  Mean: %.3f μm/s', mean(vtheta_valid))
                sprintf('  Median: %.3f μm/s', median(vtheta_valid))
                sprintf('  Std: %.3f μm/s', std(vtheta_valid))
                sprintf('  Max: %.3f μm/s', max(vtheta_valid))
                sprintf('  Min: %.3f μm/s', min(vtheta_valid))
                ''
                sprintf('Mask area: %.0f px²', areaMask(fr))
                sprintf('Centroid: (%.1f, %.1f)', xc, yc)
            };
            text(0.1, 0.9, stats_text, 'Units', 'normalized', ...
                'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'FixedWidth');
        end

        sgtitle(sprintf('Tangential Flow Analysis - Frame %d', fr), 'FontSize', 14, 'FontWeight', 'bold');

        exportgraphics(figSnap, fullfile(outDir_snapshots, sprintf('tangential_flow_analysis_fr%04d.png', fr)), ...
            'Resolution', 200);
        close(figSnap);
    end

    nSavedSnapshots = numel(1:snapshotEveryNFrames:nFrames);
    fprintf('  - Saved %d tangential flow analysis plots to %s\n', nSavedSnapshots, outDir_snapshots);
end

fprintf('All visualizations complete!\n\n');

save(fullfile(outDir,'tangential_kymo_results.mat'), ...
     'Vtheta_kymo','Npts_kymo','thetaCenters','thetaBinEdges','time_min', ...
     'centroidXY','areaMask','qcFlag','RADIUS_OF_CURVATURE', ...
     'px_per_um','dt_sec','bandOuterPx','bandInnerPx','pivMatFile','base_dir');

fprintf('Saved outputs to: %s\n', outDir);


%% ============================== FUNCTIONS =================================
function [BW, stats, poly, RADIUS, xc, yc] = make_oocyte_mask_curvature(I, ...
    sigmaBlur, threshFrac, se, polyOrder, nBoundary, theta, ...
    minAreaFrac, maxEccentric, minSolidity, centerHint)
% Curvature-based oocyte boundary reconstruction (adapted from kymograph.m).
% Methodology matches SCW_flows_curvature.m for strict physical encoding.
%
% Steps:
% 1) Coarse threshold mask (darker oocyte on brighter background)
% 2) Edge detection on coarse mask
% 3) Circle fit to get center (xc, yc) - optionally using centerHint
% 4) Polar coordinate transformation r(theta) for edge points
% 5) Polynomial fit of r(theta) with wrap-around handling (two passes)
% 6) Curvature calculation from polar curve
% 7) Reconstruct smooth boundary polygon
%
% centerHint (optional): [xc_prev, yc_prev] from previous frame for faster init

[H, W] = size(I);

% Handle optional centerHint parameter
if nargin < 11 || isempty(centerHint)
    centerHint = [W/2, H/2];  % default: image center
end

% Default outputs
BW = false(size(I));
stats = [];
poly = [];
RADIUS = nan(1, numel(theta)-1);
xc = centerHint(1); yc = centerHint(2);  % Use hint as initial guess

% 1) Coarse threshold mask
Iblur = imgaussfilt(I, sigmaBlur);
BW0 = Iblur < (threshFrac * mean2(Iblur));
BW0 = imdilate(BW0, se);
BW0 = imfill(BW0, 'holes');
BW0 = imerode(BW0, se);

% Keep largest connected component
labels = bwlabel(BW0);
if max(labels(:)) == 0
    return;  % no regions found
end
Area = zeros(1, max(labels(:)));
for i = 1:max(labels(:))
    temp = regionprops(labels == i, 'Area');
    Area(i) = temp.Area;
end
mask = labels == find(Area == max(Area), 1);
BW0 = mask;

% 2) Edge pixels of coarse mask
edges = imgradient(BW0);
edges = edges > 0;

% Extract edge pixel coordinates
xx = []; yy = [];
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        if edges(i, j) == 1
            xx = [xx j];
            yy = [yy i];
        end
    end
end

if numel(xx) < 50
    warning('Too few edge pixels for circfit; returning empty mask.');
    return;
end

% 3) Circle fit to get center (uses circfit.m from codebase)
[R, xc, yc] = circfit(xx, yy);

% 4) Compute r and angle for each edge point (polar coordinates)
nPts = numel(xx);
r = zeros(1, nPts);
angle = zeros(1, nPts);

for kk = 1:nPts
    r(kk) = norm([xx(kk) - xc, yy(kk) - yc]);

    % Angle calculation matching kymograph.m methodology
    if yy(kk) > yc
        angle(kk) = acos(dot(([xx(kk) - xc, yy(kk) - yc]) / norm([xx(kk) - xc, yy(kk) - yc]), [1 0]));
    else
        angle(kk) = 2*pi - acos(dot(([xx(kk) - xc, yy(kk) - yc]) / norm([xx(kk) - xc, yy(kk) - yc]), [1 0]));
    end
end

% 5) FIRST PASS: fit r(angle) and fill middle half of RR
RRrow = nan(1, nBoundary);

param = polyfit(angle, r, polyOrder);
xgrid = linspace(0, 2*pi, nBoundary);
y1 = polyval(param, xgrid);

r_theta_p = polyder(param);
r_2theta_p = polyder(r_theta_p);
r_theta = polyval(r_theta_p, xgrid);
r_2theta = polyval(r_2theta_p, xgrid);

% Curvature calculation for middle half (bins 26-75 for 101 theta bins)
for ii = (length(theta)-1)/2 - (length(theta)-1)/4 : (length(theta)-1)/2 + (length(theta)-1)/4
    sel = xgrid > theta(ii) & xgrid < theta(ii+1);
    if any(sel)
        num = ((y1(sel).^2 + r_theta(sel).^2).^(3/2));
        den = abs(y1(sel).^2 + 2*r_theta(sel).^2 - y1(sel).*r_2theta(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

RRrow(126:375) = y1(126:375);

% 6) SECOND PASS: sort and shift by pi to handle wrap-around
[angle2, Iord] = sort(angle);
r2 = r(Iord);

angle2 = angle2 + pi;
angle2(angle2 > 2*pi) = angle2(angle2 > 2*pi) - 2*pi;

param2 = polyfit(angle2, r2, polyOrder);
x2 = linspace(0, 2*pi, nBoundary);
y2 = polyval(param2, x2);

r_theta_p2 = polyder(param2);
r_2theta_p2 = polyder(r_theta_p2);
r_theta2 = polyval(r_theta_p2, x2);
r_2theta2 = polyval(r_2theta_p2, x2);

x2 = x2 - pi;
x2(x2 < 0) = 2*pi + x2(x2 < 0);

y2 = circshift(y2, nBoundary/2);

RRrow(1:125) = y2(1:125);
RRrow(376:500) = y2(376:500);

% Curvature for remaining quarters
for ii = 1:(length(theta)-1)/2 - (length(theta)-1)/4
    sel = x2 > theta(ii) & x2 < theta(ii+1);
    if any(sel)
        num = ((y2(sel).^2 + r_theta2(sel).^2).^(3/2));
        den = abs(y2(sel).^2 + 2*r_theta2(sel).^2 - y2(sel).*r_2theta2(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

for ii = (length(theta)-1)/2 + (length(theta)-1)/4 : length(theta)-1
    sel = x2 > theta(ii) & x2 < theta(ii+1);
    if any(sel)
        num = ((y2(sel).^2 + r_theta2(sel).^2).^(3/2));
        den = abs(y2(sel).^2 + 2*r_theta2(sel).^2 - y2(sel).*r_2theta2(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

% 7) Reconstruct final smooth boundary from RR
xfinal = linspace(0, 2*pi, nBoundary);
XX = RRrow .* cos(xfinal) + xc;
YY = RRrow .* sin(xfinal) + yc;

% Clamp to image bounds for poly2mask stability
XX = min(max(XX, 1), W);
YY = min(max(YY, 1), H);

% Build final mask from reconstructed smooth boundary
BW = poly2mask(XX, YY, H, W);
BW = imfill(BW, 'holes');

poly = [XX(:) YY(:)];

% QC: check plausibility
st = regionprops(BW, 'Area', 'Eccentricity', 'Centroid', 'Solidity');
if isempty(st)
    BW = false(size(BW));
    stats = [];
    poly = [];
    return;
end
[~, ii] = max([st.Area]);
stats = st(ii);

% QC thresholds
imgArea = numel(I);
if stats.Area < minAreaFrac*imgArea || stats.Eccentricity > maxEccentric || stats.Solidity < minSolidity
    BW = false(size(BW));
    stats = [];
    poly = [];
    return;
end

end

function [vBins, nCount, dbg] = tangential_from_boundary_curvature( ...
    X, Y, U, V, BW, poly, centroid, bandOuterPx, bandInnerPx, nThetaBins, thetaBinEdges, ...
    px_per_um, dt_sec, pivVelUnit, nDenseSpline)
% Calculate tangential velocity field using curvature-based boundary.
% STRICT PHYSICAL ENCODING:
% - Tangent vectors computed from smooth polynomial boundary
% - Velocities converted to physical units (um/s)
% - Tangential component: v_tangent = v · t̂ (dot product with unit tangent)
% - Binned by angular position theta around oocyte centroid

% Convert velocities to um/s (strict physical units)
[U_um_s, V_um_s] = convertVelToUmPerSec(U, V, pivVelUnit, px_per_um, dt_sec);

% Use polygon boundary from curvature reconstruction
if isempty(poly)
    vBins = nan(1,nThetaBins);
    nCount = zeros(1,nThetaBins);
    dbg = struct('xDense',[],'yDense',[],'sampleX',[],'sampleY',[]);
    return;
end

xb = poly(:,1);
yb = poly(:,2);

% Close boundary explicitly
if ~isequal([xb(1) yb(1)], [xb(end) yb(end)])
    xb = [xb; xb(1)];
    yb = [yb; yb(1)];
end

% Arclength parameter s
ds = hypot(diff(xb), diff(yb));
s  = [0; cumsum(ds)];
L  = s(end);

if L < 50
    vBins = nan(1,nThetaBins);
    nCount = zeros(1,nThetaBins);
    dbg = struct('xDense',xb,'yDense',yb,'sampleX',[],'sampleY',[]);
    return;
end

% Simple finite difference tangent calculation (centered differences)
% Faster and doesn't require Curve Fitting Toolbox
% Densely sample boundary using linear interpolation
tSample = linspace(0, 1, nDenseSpline);
xDense = interp1(linspace(0,1,numel(xb)), xb, tSample, 'linear');
yDense = interp1(linspace(0,1,numel(yb)), yb, tSample, 'linear');

nPts = numel(xDense);
tx = zeros(nPts, 1);
ty = zeros(nPts, 1);

% Centered finite difference with periodic boundary conditions
for i = 1:nPts
    i_prev = mod(i-2, nPts) + 1;  % wrap around for i=1
    i_next = mod(i, nPts) + 1;     % wrap around for i=nPts

    % Tangent = (next_point - prev_point) / 2
    % This approximates dr/ds using centered difference
    tx(i) = (xDense(i_next) - xDense(i_prev)) / 2;
    ty(i) = (yDense(i_next) - yDense(i_prev)) / 2;
end

% Normalize tangent to unit vector (strict physical encoding: t̂)
tN = hypot(tx,ty) + eps;
tx = tx./tN;
ty = ty./tN;

% Define cortical band using distance-to-perimeter inside BW
% This ensures we only sample velocities within the cortical region
per = bwperim(BW);
D = bwdist(per);

% Clamp indices for D lookup
Xi = clamp(round(X), 1, size(D,2));
Yi = clamp(round(Y), 1, size(D,1));
lin = sub2ind(size(D), Yi, Xi);

inside = BW(lin) & isfinite(U_um_s) & isfinite(V_um_s);
inBand = inside & (D(lin) >= bandOuterPx) & (D(lin) <= bandInnerPx);

if ~any(inBand(:))
    vBins = nan(1,nThetaBins);
    nCount = zeros(1,nThetaBins);
    dbg = struct('xDense',xDense,'yDense',yDense,'sampleX',[],'sampleY',[]);
    return;
end

xq = X(inBand); yq = Y(inBand);
uq = U_um_s(inBand); vq = V_um_s(inBand);

% Nearest point on dense boundary samples (brute-force; PIV grids are small)
idxNearest = nearest_dense_points(xq, yq, xDense, yDense);

tqx = tx(idxNearest);
tqy = ty(idxNearest);

% STRICT PHYSICAL ENCODING: tangential component = v · t̂
% where v = (u, v) is velocity vector, t̂ = (tx, ty) is unit tangent
vT = uq.*tqx + vq.*tqy;  % tangential velocity [um/s]

% Theta coordinate assigned from nearest boundary point relative to centroid
% This maintains angular position consistency with curvature calculation
cx = centroid(1); cy = centroid(2);
th = atan2(yDense(idxNearest) - cy, xDense(idxNearest) - cx);
th = wrapTo2Pi(th);

% Bin into theta bins (angular averaging)
vBins = nan(1,nThetaBins);
nCount = zeros(1,nThetaBins);
binIdx = discretize(th, thetaBinEdges);
for k=1:nThetaBins
    m = (binIdx == k);
    nCount(k) = sum(m);
    if any(m), vBins(k) = mean(vT(m), 'omitnan'); end
end

dbg = struct();
dbg.xDense = xDense; dbg.yDense = yDense;
dbg.sampleX = xq; dbg.sampleY = yq;

end

function idx = nearest_dense_points(xq, yq, xd, yd)
idx = zeros(size(xq));
for i=1:numel(xq)
    dx = xd - xq(i);
    dy = yd - yq(i);
    [~, idx(i)] = min(dx.*dx + dy.*dy);
end
end

function v = clamp(v, lo, hi)
v = max(lo, min(hi, v));
end

function [U_um_s, V_um_s] = convertVelToUmPerSec(U, V, pivVelUnit, px_per_um, dt_sec)
switch lower(strtrim(pivVelUnit))
    case 'px_per_frame'
        U_um_s = (U / px_per_um) / dt_sec;
        V_um_s = (V / px_per_um) / dt_sec;
    case 'px_per_sec'
        U_um_s = (U / px_per_um);
        V_um_s = (V / px_per_um);
    otherwise
        error('Unknown pivVelUnit: %s', pivVelUnit);
end
end

function [Xc, Yc, Uc, Vc] = pickPIVFields(S)
% Try common PIVlab exports
if isfield(S,'X') && isfield(S,'Y')
    Xc = S.X; Yc = S.Y;
elseif isfield(S,'x') && isfield(S,'y')
    Xc = S.x; Yc = S.y;
else
    error('Cannot find X/Y in PIV .mat (expected X,Y or x,y).');
end

if isfield(S,'U') && isfield(S,'V')
    Uc = S.U; Vc = S.V;
elseif isfield(S,'u') && isfield(S,'v')
    Uc = S.u; Vc = S.v;
elseif isfield(S,'u_original') && isfield(S,'v_original')
    Uc = S.u_original; Vc = S.v_original;
else
    error('Cannot find U/V in PIV .mat (expected U,V or u,v or u_original,v_original).');
end

if ~iscell(Uc), Uc={Uc}; end
if ~iscell(Vc), Vc={Vc}; end
if ~iscell(Xc), Xc={Xc}; end
if ~iscell(Yc), Yc={Yc}; end

% replicate X/Y if single but U/V multi-frame
if numel(Xc)==1 && numel(Uc)>1, Xc = repmat(Xc, size(Uc)); end
if numel(Yc)==1 && numel(Uc)>1, Yc = repmat(Yc, size(Uc)); end
end

function cmap = redblue(n)
% REDBLUE  Diverging red-white-blue colormap for signed data
% Blue (negative) -> White (zero) -> Red (positive)
if nargin < 1, n = 256; end

r = linspace(0, 1, n)';
cmap = [r, 1-abs(2*r-1), 1-r];  % Simple red-white-blue

% Alternative: use MATLAB's built-in if available
if exist('brewermap', 'file')
    cmap = brewermap(n, 'RdBu');
    cmap = flipud(cmap);  % flip so red=positive
end
end